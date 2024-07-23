__author__ = 'chris'

import subprocess
import os
import re

'''
This is a pre processor class for the Tuxedo Suite.
'''

__all__ = [
    "TPreProcessor",
]

class TPreProcessor():

    VERY_SENSITIVE = 4
    SENSITIVE = 3
    FAST = 2
    VERY_FAST = 1

    def __init__(self, refgenome_name, rawdata_filepath, transcripts_name, transcripts_filepath=None,
                 rawdata_paired2_filepath=None, refgenome_filepath=None, data_outdir='data/out/'):

        """
        :param refgenome_name: String
        :param rawdata_filepath: String
        :param transcripts_name: String
        :param rawdata_paired2_filepath: String
        :param transcripts_filepath: String
        :param refgenome_filepath: String
        :param data_outdir: String
        :return: None
        """
        self._dataset = []
        self._refgenome_name = refgenome_name
        self._data_outdir = data_outdir
        self._rawdata_filepath = rawdata_filepath
        self._transcripts_name = transcripts_name

        if (rawdata_paired2_filepath is not None):
            self._isPairedEnd = True
            self._rawdata_paired2_filepath = rawdata_paired2_filepath

        if (refgenome_filepath is not None):
            self._refgenome_filepath = refgenome_filepath
            self.indexRefGenome(self._refgenome_filepath)

        if (transcripts_filepath is not None):
            self._transcripts_filepath = transcripts_filepath
            self.indexRefTranscript(self._transcripts_filepath, transcripts_name)


    def indexRefGenome(self, refgenome_filepath):
        """Index the reference genome using Bowtie Build
        :param refgenome_filepath: String
        """

        subprocess.call(["bowtie2-build", refgenome_filepath, self._refgenome_name])


    def indexRefTranscript(self, reftranscripts_filepath):
        """ Index the reference transcripts using Bowtie Build
        :param reftranscripts_filepath: String
        """

        subprocess.call(["bowtie2-build", reftranscripts_filepath, self._transcripts_name])


    def alignSingleReadsToGenome(self, output_filename):
        """ Align the given single read fastq files to the reference genome using Bowtie2
        :param output_filename: String
        """
        if (not self._isPairedEnd):
            subprocess.call(
                ["bowtie2", self._refgenome_name, self._rawdata_filepath, "-S", self._data_outdir + output_filename])
        else:
            raise Exception("Class Instance is set with paried end reads. You can not align single reads.")


    def alignPairedReadsToGenome(self, output_filename):
        """Align the given paired end reads to the reference genome using Bowtie2
        :param output_filename: String
        """
        if (self._isPairedEnd):
            subprocess.call(["bowtie2", "-x", self._refgenome_name, "-1", self._rawdata_filepath, "-2",
                             self._rawdata_paired2_filepath, "-S", self._data_outdir + output_filename])
        else:
            raise Exception("Class Instance is set with single end reads. You can not align paried reads.")


    def alignSingleReadsToTrans(self, output_filename):
        """Align the given single end reads to the reference transcripts using Bowtie2
        :param output_filename: String
        """
        if (not self._isPairedEnd):
            subprocess.call(
                ["bowtie2", self._transcripts_name, self._rawdata_filepath, "-S", self._data_outdir + output_filename])
        else:
            raise Exception("Class Instance is set with paried end reads. You can not align single reads.")


    def alignPairedReadsToTrans(self, output_filename):
        """Align the given paired end reads to the reference transcripts using Bowtie2
        :param output_filename: String
        """

        if (self._isPairedEnd):
            subprocess.call(["bowtie2", "-x", self._transcripts_name, "-1", self._rawdata_filepath, "-2",
                             self._rawdata_paired2_filepath, "-S", self._data_outdir + output_filename])
        else:
            raise Exception("Class Instance is set with single end reads. You can not align paried reads.")


    def samToBam(self, sam_filepath, output_filename, remove_sam=False):
        """Convert the sam file produced by alignment into a bam file
        :param remove_sam: Boolean
        :param output_filename: String
        """
        subprocess.call(["samtools", "view", "-bSo", self._data_outdir + output_filename, self._data_outdir + sam_filepath])
        if (remove_sam):
            subprocess.call(["rm", sam_filepath])

    def sortBam(self, bam_filepath, output_filename):
        """Sort the bam file by chromosomal coordinate
        :param bam_filepath: String
        :param output_filename: String
        """
        subprocess.call(["samtools", "sort", self._data_outdir + bam_filepath, self._data_outdir + output_filename])


    def indexBam(self, sbam_filepath):
        """Index the bam file
        :param sbam_filepath: String
        """
        subprocess.call(["samtools", "index", self._data_outdir + sbam_filepath])


    def __getmetrics(self, sorted_bam, out_text):
        """get the mean_insert_size and standard_deviation to give to the tophat command
        :param sorted_bam: String
        :param out_text: String
        """
        histogram_file = "test_hist.pdf"
        try:
            subprocess.call(
                ["picard-tools", "CollectInsertSizeMetrics", "INPUT=", self._data_outdir + sorted_bam, "HISTOGRAM_FILE=", self._data_outdir + histogram_file, "OUTPUT=", self._data_outdir + out_text])
            stdin = os.popen("head -8 " + self._data_outdir + out_text + " | tail -2 | cut -f 5,6")
            metrics = stdin.read()
            datalist = re.split('\n|\t', metrics)
            datalist = datalist[:-1]
            mein = float(datalist[2])
            stde = float(datalist[3])
            mein = round(mein)
            stde = round(stde)
            mein = int(mein)
            stde = int(stde)

            if (self._isPairedEnd):
                sin = os.popen("head -2 " + self._rawdata_filepath + " | tail -1 | tr -d '\n' | wc -c")
                readLength = sin.read()
                readLength = int(readLength)
                mein -= 2 * readLength  #Subtract twice the readlength from the mean_insert size

            self._dataset = [mein, stde]
        except Exception as e:
            print(e)



    def mapReadsWithTopHat(self, sorted_bam, out_text, sensetivity, output_directory):
        """Map the reads using tophat
        :param sorted_bam: String
        :param out_text: String
        :param sensetivity: Int
        :param output_directory: String
        """

        self.__getmetrics(sorted_bam, out_text)

        s_command = ""
        if (sensetivity < self.VERY_FAST or sensetivity > self.VERY_SENSITIVE):
            raise Exception("Sensetivity must have a value >=1 and <=4")

        if (sensetivity == self.VERY_FAST):
            #very fast
            s_command = "--b2-very-fast"
        elif (sensetivity == self.FAST):
            #fast
            s_command = "--b2-fast"
        elif (sensetivity == self.SENSITIVE):
            #sensitive
            s_command = "--b2-sensitive"
        elif (sensetivity == self.VERY_SENSITIVE):
            #very sensitive
            s_command = "--b2-very-sensitive"

        if (self._isPairedEnd):
            mein = str(self._dataset[0])
            stde = str(self._dataset[1])
            subprocess.call(["tophat", "-o", self._data_outdir + output_directory, s_command, "-r", mein, "--mate-std-dev", stde, self._refgenome_name, self._rawdata_filepath, self._rawdata_paired2_filepath])
        else:
            raise Exception("This tophat method can only be used for paired end data")


    def assembleTranscripts(self, tophat_outdir, output_directory="cufflinks"):
        """  Assemble the transcripts using Cufflinks
        :param tophat_outdir: String
        :param output_directory: String
        """
        subprocess.call(["cufflinks", "-o", self._data_outdir + output_directory, self._data_outdir + tophat_outdir + "/accepted_hits.bam"])

#have a switch on all methods for users to overide default location for finding data.

    def mergeTranscripts(self, tfile_list, reference_annotation=None, output_directory="cuffmerge"):
        """Merge the transcripts using cufffmerge
        :param tfile_list: List of Strings
        :param output_directory: String
        :param reference_annotation: String
        """
        tran_path = self._data_outdir + "transcripts.txt"

        subprocess.call(["touch", tran_path])
        f = open(tran_path, "w")
        for path in tfile_list:
            f.write(self._data_outdir + path + "/transcripts.gtf\n")

        f.close()

        if (reference_annotation is not None):
            subprocess.call(["cuffmerge", "-o", self._data_outdir + output_directory, "-g", reference_annotation, tran_path])
        else:
            subprocess.call(["cuffmerge", "-o", self._data_outdir + output_directory, tran_path])
        subprocess.call(["rm", tran_path])


    def cuffdiff(self, transcripts_filepath, bam_files, output_directory="cuffdiff"):
        """
        Perform differental analysis on the two (or more) bam files.
        :param transcripts_filepath: String
        :param bam_files: list of Strings
        :param output_directory: String
        """

        command = "cuffdiff " + "-o " + self._data_outdir + output_directory + " " + self._data_outdir + transcripts_filepath

        for path in bam_files:

            command = command + " " + self._data_outdir + path

        command = command.split(" ")
        subprocess.call(command)




