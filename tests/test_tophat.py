__author__ = 'chris'


from tophat_package.tuxedpreprocessing import TPreProcessor

#tpp = TPreProcessor("data/sequencing/rna_seq/hg19", "data/sequencing/rna_seq/brain_1.fastq", "data/sequencing/rna_seq/hg19_transcripts",
                                #rawdata_paired2_filepath="data/sequencing/rna_seq/brain_2.fastq", data_outdir="data/sequencing/results/")

#tpp.alignPairedReadsToTrans("brain.sam")
#tpp.samToBam("brain.sam", "brain.bam", True)
#tpp.sortBam("brain.bam", "brain.sorted")
#tpp.mapReadsWithTopHat("brain.sorted.bam", "brain_out.txt", TPreProcessor.VERY_SENSITIVE, "brain_tophat")
#tpp.assembleTranscripts("brain_tophat", "brain_cufflinks")


tpp = TPreProcessor("data/sequencing/rna_seq/hg19", "data/sequencing/rna_seq/adrenal_1.fastq", "data/sequencing/rna_seq/hg19_transcripts",
                                rawdata_paired2_filepath="data/sequencing/rna_seq/adrenal_2.fastq", data_outdir="data/sequencing/results/")

#tpp.alignPairedReadsToTrans("adrenal.sam")
#tpp.samToBam("adrenal.sam", "adrenal.bam", True)
#tpp.sortBam("adrenal.bam", "adrenal.sorted")
#tpp.mapReadsWithTopHat("adrenal.sorted.bam", "adrenal_out.txt", TPreProcessor.VERY_SENSITIVE, "adrenal_tophat")
#tpp.assembleTranscripts("adrenal_tophat", "adrenal_cufflinks")

#tpp.mergeTranscripts(["brain_cufflinks", "adrenal_cufflinks"], "data/sequencing/rna_seq/hg19_chr19_gene_annotation.gtf", output_directory="cuffmerge")
tpp.cuffdiff("cuffmerge/merged.gtf", ["adrenal_tophat/accepted_hits.bam", "brain_tophat/accepted_hits.bam"], "cuffdiff")