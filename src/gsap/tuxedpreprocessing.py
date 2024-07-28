"""
This module implements functions for analysing RNA-Seq data. It is used for pre-processing raw single-end short-read sequence data initially in FASTQ format into other data forms as required by tools such as samtools, picard-tools, and the trio of bowtie2, cufflinks, and tophat (so called the tuxedo suite).

Author: Chris Talor, Sunny Park

"""

import os
import re
import subprocess

__all__ = [
    "TPreProcessor",
]


class TPreProcessor:
    VERY_SENSITIVE = 4
    SENSITIVE = 3
    FAST = 2
    VERY_FAST = 1

    def __init__(
        self,
        refgenome_name: str,
        rawdata_filepath: str,
        transcripts_name: str,
        transcripts_filepath: str | None = None,
        rawdata_paired2_filepath: str | None = None,
        refgenome_filepath: str | None = None,
        data_outdir: str = "data/out/",
    ) -> None:
        """Constructor method."""
        self._dataset: list[int] = []
        self._refgenome_name = refgenome_name
        self._data_outdir = data_outdir
        self._rawdata_filepath = rawdata_filepath
        self._transcripts_name = transcripts_name

        if rawdata_paired2_filepath is not None:
            self._isPairedEnd = True
            self._rawdata_paired2_filepath = rawdata_paired2_filepath

        if refgenome_filepath is not None:
            self._refgenome_filepath = refgenome_filepath
            self.indexRefGenome(self._refgenome_filepath)

        if transcripts_filepath is not None:
            self._transcripts_filepath = transcripts_filepath
            self.indexRefTranscript(self._transcripts_filepath)

    def indexRefGenome(self, refgenome_filepath: str) -> None:
        """Index the reference genome using Bowtie Build."""
        subprocess.call(
            ["/usr/bin/bowtie2-build", refgenome_filepath, self._refgenome_name]
        )

    def indexRefTranscript(self, reftranscripts_filepath: str) -> None:
        """Index the reference transcripts using Bowtie Build."""
        subprocess.call(
            ["/usr/bin/bowtie2-build", reftranscripts_filepath, self._transcripts_name]
        )

    def alignSingleReadsToGenome(self, output_filename: str) -> None:
        """Align the given single read fastq files to the reference genome using Bowtie2."""
        if not self._isPairedEnd:
            subprocess.call(
                [
                    "/usr/bin/bowtie2",
                    self._refgenome_name,
                    self._rawdata_filepath,
                    "-S",
                    self._data_outdir + output_filename,
                ]
            )
        else:
            raise Exception(
                "Class Instance is set with paried end reads. You can not align single reads."
            )

    def alignPairedReadsToGenome(self, output_filename: str) -> None:
        """Align the given paired end reads to the reference genome using Bowtie2."""
        if self._isPairedEnd:
            subprocess.call(
                [
                    "/usr/bin/bowtie2",
                    "-x",
                    self._refgenome_name,
                    "-1",
                    self._rawdata_filepath,
                    "-2",
                    self._rawdata_paired2_filepath,
                    "-S",
                    self._data_outdir + output_filename,
                ]
            )
        else:
            raise Exception(
                "Class Instance is set with single end reads. You can not align paried reads."
            )

    def alignSingleReadsToTrans(self, output_filename: str) -> None:
        """Align the given single end reads to the reference transcripts using Bowtie2."""
        if not self._isPairedEnd:
            subprocess.call(
                [
                    "/usr/bin/bowtie2",
                    self._transcripts_name,
                    self._rawdata_filepath,
                    "-S",
                    self._data_outdir + output_filename,
                ]
            )
        else:
            raise Exception(
                "Class Instance is set with paried end reads. You can not align single reads."
            )

    def alignPairedReadsToTrans(self, output_filename: str) -> None:
        """Align the given paired end reads to the reference transcripts using Bowtie2."""
        if self._isPairedEnd:
            subprocess.call(
                [
                    "/usr/bin/bowtie2",
                    "-x",
                    self._transcripts_name,
                    "-1",
                    self._rawdata_filepath,
                    "-2",
                    self._rawdata_paired2_filepath,
                    "-S",
                    self._data_outdir + output_filename,
                ]
            )
        else:
            raise Exception(
                "Class Instance is set with single end reads. You can not align paried reads."
            )

    def samToBam(
        self, sam_filepath: str, output_filename: str, remove_sam: bool = False
    ) -> None:
        """Convert the sam file produced by alignment into a bam file."""
        subprocess.call(
            [
                "/usr/bin/samtools",
                "view",
                "-bSo",
                self._data_outdir + output_filename,
                self._data_outdir + sam_filepath,
            ]
        )
        if remove_sam:
            subprocess.call(["/usr/bin/rm", sam_filepath])

    def sortBam(self, bam_filepath: str, output_filename: str) -> None:
        """Sort the bam file by chromosomal coordinate."""
        subprocess.call(
            [
                "/usr/bin/samtools",
                "sort",
                self._data_outdir + bam_filepath,
                self._data_outdir + output_filename,
            ]
        )

    def indexBam(self, sbam_filepath: str) -> None:
        """Index the bam file."""
        subprocess.call(
            ["/usr/bin/samtools", "index", self._data_outdir + sbam_filepath]
        )

    def __getmetrics(self, sorted_bam: str, out_text: str) -> None:
        """Get the mean_insert_size and standard_deviation to give to the tophat command."""
        histogram_file = "test_hist.pdf"
        try:
            subprocess.call(
                [
                    "/usr/bin/picard-tools",
                    "CollectInsertSizeMetrics",
                    "INPUT=",
                    self._data_outdir + sorted_bam,
                    "HISTOGRAM_FILE=",
                    self._data_outdir + histogram_file,
                    "OUTPUT=",
                    self._data_outdir + out_text,
                ]
            )
            stdin = os.popen(
                "head -8 " + self._data_outdir + out_text + " | tail -2 | cut -f 5,6"
            )
            metrics = stdin.read()
            datalist = re.split("\n|\t", metrics)
            datalist = datalist[:-1]
            mein = float(datalist[2])
            stde = float(datalist[3])
            mein = round(mein)
            stde = round(stde)
            mein = int(mein)
            stde = int(stde)

            if self._isPairedEnd:
                sin = os.popen(
                    "head -2 "
                    + self._rawdata_filepath
                    + " | tail -1 | tr -d '\n' | wc -c"
                )
                readLengthStr = sin.read()
                readLength = int(readLengthStr)
                mein -= (
                    2 * readLength
                )  # Subtract twice the readlength from the mean_insert size

            self._dataset = [mein, stde]
        except Exception as e:
            print(e)

    def mapReadsWithTopHat(
        self, sorted_bam: str, out_text: str, sensitivity: int, output_directory: str
    ) -> None:
        """Map the reads using tophat."""
        self.__getmetrics(sorted_bam, out_text)

        s_command = ""
        if sensitivity < self.VERY_FAST or sensitivity > self.VERY_SENSITIVE:
            raise Exception("Sensitivity must have a value >=1 and <=4")

        if sensitivity == self.VERY_FAST:
            # very fast
            s_command = "--b2-very-fast"
        elif sensitivity == self.FAST:
            # fast
            s_command = "--b2-fast"
        elif sensitivity == self.SENSITIVE:
            # sensitive
            s_command = "--b2-sensitive"
        elif sensitivity == self.VERY_SENSITIVE:
            # very sensitive
            s_command = "--b2-very-sensitive"

        if self._isPairedEnd:
            mein = str(self._dataset[0])
            stde = str(self._dataset[1])
            subprocess.call(
                [
                    "/usr/bin/tophat",
                    "-o",
                    self._data_outdir + output_directory,
                    s_command,
                    "-r",
                    mein,
                    "--mate-std-dev",
                    stde,
                    self._refgenome_name,
                    self._rawdata_filepath,
                    self._rawdata_paired2_filepath,
                ]
            )
        else:
            raise Exception("This tophat method can only be used for paired end data")

    def assembleTranscripts(
        self, tophat_outdir: str, output_directory: str = "cufflinks"
    ) -> None:
        """Assemble the transcripts using Cufflinks."""
        subprocess.call(
            [
                "/usr/bin/cufflinks",
                "-o",
                self._data_outdir + output_directory,
                self._data_outdir + tophat_outdir + "/accepted_hits.bam",
            ]
        )

    # have a switch on all methods for users to override default location
    # for finding data.

    def mergeTranscripts(
        self,
        tfile_list: list[str],
        reference_annotation: str | None,
        output_directory: str = "cuffmerge",
    ) -> None:
        """Merge the transcripts using cufffmerge."""
        tran_path = self._data_outdir + "transcripts.txt"

        subprocess.call(["/usr/bin/touch", tran_path])
        f = open(tran_path, "w")
        for path in tfile_list:
            f.write(self._data_outdir + path + "/transcripts.gtf\n")

        f.close()

        if reference_annotation is not None:
            subprocess.call(
                [
                    "/usr/bin/cuffmerge",
                    "-o",
                    self._data_outdir + output_directory,
                    "-g",
                    reference_annotation,
                    tran_path,
                ]
            )
        else:
            subprocess.call(
                [
                    "/usr/bin/cuffmerge",
                    "-o",
                    self._data_outdir + output_directory,
                    tran_path,
                ]
            )
        subprocess.call(["/usr/bin/rm", tran_path])

    def cuffdiff(
        self,
        transcripts_filepath: str,
        bam_files: list[str],
        output_directory: str = "cuffdiff",
    ) -> None:
        """Perform differental analysis on the two (or more) bam files."""
        command = (
            "cuffdiff "
            + "-o "
            + self._data_outdir
            + output_directory
            + " "
            + self._data_outdir
            + transcripts_filepath
        )

        for path in bam_files:
            command = command + " " + self._data_outdir + path

        commandstr = command.split(" ")
        subprocess.call(commandstr)
