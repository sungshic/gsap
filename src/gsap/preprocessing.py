__author__ = 'spark'

import os
import subprocess

__all__ = [
    "PreProcessor",
]

class PreProcessor():
    def __init__(self, home_dir, refgenome_name, rawdata_filepath, rawdata_paired2_filepath=None, refgenome_filepath=None, data_outdir='data/out/', tmpdir='data/results/tmp'):
        self._rawdata_filepath = rawdata_filepath
        self._home_dir = home_dir
        if (rawdata_paired2_filepath != None):
            self._isPairedEnd = True
            self._rawdata_paired2_filepath = rawdata_paired2_filepath
        if (refgenome_filepath != None):
            self._refgenome_filepath = refgenome_filepath
            self.indexRefGenome(self._refgenome_filepath, refgenome_name)
        self._refgenome_name = refgenome_name
        self._data_outdir = data_outdir
        self._tmpdir = tmpdir

    def qualityCheck(self, fastq_filepath):
        subprocess.call([self._home_dir+'toolset/FastQC/fastqc', fastq_filepath])

    def trimSeqReads(self):
        fastq_filepath_for = self._rawdata_filepath
        fastq_filepath_rev = self._rawdata_paired2_filepath
        out_forward_paired = self._data_outdir + 'out_for_paired.fq.gz'
        out_forward_unpaired = self._data_outdir + 'out_for_unpaired.fq.gz'
        out_rev_paired = self._data_outdir + 'out_rev_paired.fq.gz'
        out_rev_unpaired = self._data_outdir + 'out_rev_unpaired.fq.gz'

        try:
            subprocess.call(['java', '-jar', self._home_dir+'toolset/Trimmomatic-0.32/trimmomatic-0.32.jar', 'PE', '-phred33', fastq_filepath_for, fastq_filepath_rev, out_forward_paired, out_forward_unpaired, out_rev_paired, out_rev_unpaired, 'LEADING:28', 'TRAILING:28', 'HEADCROP:15', 'MINLEN:36']) #, sort_order])
            # replace the raw files with trimmed sequences
            self._rawdata_filepath = out_forward_paired
            self._rawdata_paired2_filepath = out_rev_paired
        except Exception as e:
            print(e)

    def assembleDeNovo(self): #, kmer_range):
        out_forward_paired = self._data_outdir + 'out_for_paired.fq.gz'
        out_rev_paired = self._data_outdir + 'out_rev_paired.fq.gz'
        subprocess.call(['python', self._home_dir+'toolset/SPAdes-4.0.0-Linux/bin/spades.py', '--pe1-1', out_forward_paired, '--pe1-2', out_rev_paired, '-o', self._data_outdir+ 'denovo/spades/output'])
    
    def indexRefGenome(self, refgenome_filepath, refgenome_name):
        subprocess.call([self._home_dir+'toolset/bowtie2/bowtie2-build', refgenome_filepath, refgenome_name])

    def alignPairedEndReads(self, output_filename):
        self._alignPairedEndReads(output_filename, self._refgenome_name, self._rawdata_filepath, self._rawdata_paired2_filepath)

    def _alignPairedEndReads(self, output_filename, ref_name, reads1_filepath, reads2_filepath):
        subprocess.call([self._home_dir+'toolset/bowtie2/bowtie2', '-x', ref_name, '-1', reads1_filepath, '-2', reads2_filepath, '-S', self._data_outdir+output_filename])

    def sortSAMIntoBAM(self, insam_filepath, outbam_filename, sort_order='coordiante'):
        try:
            subprocess.call(['java', '-Xmx4g', '-Djava.io.tmpdir='+self._tmpdir, '-jar', self._home_dir+'toolset/picard-tools/SortSam.jar', 'INPUT=', insam_filepath, 'OUTPUT=', self._data_outdir+outbam_filename, 'SORT_ORDER=coordinate']) #, sort_order])
        except Exception as e:
            print(e)

    def sortBAMIntoBAM(self, inbam_filepath, outbam_filename, sort_order='coordiante'):
        try:
            subprocess.call(['java', '-Xmx2g', '-Djava.io.tmpdir='+self._tmpdir, '-jar', self._home_dir+'toolset/picard-tools/SortSam.jar', 'INPUT=', inbam_filepath, 'OUTPUT=', self._data_outdir+outbam_filename, 'SORT_ORDER=coordinate']) #, sort_order])
        except Exception as e:
            print(e)



    def markDuplicates(self, inbam_filepath, outbam_filename, outmetrics_filename):
        try:
            subprocess.call(['java', '-Xmx2g', '-Djava.io.tmpdir='+self._tmpdir, '-jar', self._home_dir+'toolset/picard-tools/MarkDuplicates.jar', 'INPUT=', inbam_filepath, 'OUTPUT=', self._data_outdir+outbam_filename, 'METRICS_FILE=', self._data_outdir+outmetrics_filename])
        except Exception as e:
            print(e)

    def addReadGroup(self, inbam_filepath, outbam_filename):
        try:
            subprocess.call(['java', '-Xmx2g', '-Djava.io.tmpdir='+self._tmpdir, 
                             '-jar', self._home_dir+'toolset/picard-tools/AddOrReplaceReadGroups.jar',
                             'INPUT=', inbam_filepath,
                             'OUTPUT=', self._data_outdir+outbam_filename,
                             'RGID=', 'group1',
                             'RGLB=', 'lib1',
                             'RGPL=', 'illumina',
                             'RGPU=', 'unit1',
                             'RGSM=', 'sample1'])
        except Exception as e:
            print(e)


    def buildBamIndex(self, inbam_filepath):
        try:
            subprocess.call(['java', '-Xmx2g', '-Djava.io.tmpdir='+self._tmpdir, '-jar', self._home_dir+'toolset/picard-tools/BuildBamIndex.jar', 'INPUT=', inbam_filepath])
        except Exception as e:
            print(e)

    def buildRefFastaIndex(self, inreffasta_filepath, outdict_filepath):
        try:
            subprocess.call(['rm', outdict_filepath])
            # create a dictionary file from fasta
            subprocess.call(['java', '-Xmx2g', '-Djava.io.tmpdir='+self._tmpdir, '-jar', self._home_dir+'toolset/picard-tools/CreateSequenceDictionary.jar',
                             'R=', inreffasta_filepath,
                             'O=', outdict_filepath])
            # create a fasta index file
            subprocess.call([self._home_dir+'toolset/samtools-1.2/samtools', 'faidx', inreffasta_filepath])
        except Exception as e:
            print(e)

    def removeUnmappedReads(self, inbam_filepath, mapped_outbam_filepath, unmapped_outbam_filepath):
        try:
            with open(mapped_outbam_filepath, 'w') as f_out:
                subprocess.call([self._home_dir+'toolset/samtools-1.2/samtools', 'view', '-hbF', '4', inbam_filepath], stdout=f_out)
            with open(unmapped_outbam_filepath, 'w') as f_out:
                subprocess.call([self._home_dir+'toolset/samtools-1.2/samtools', 'view', '-hbf', '4', inbam_filepath], stdout=f_out)
        except Exception as e:
            print(e)

