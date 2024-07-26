__author__ = 'spark'

import os
import sh
import subprocess
from gsap import PreProcessor, VariantCaller, SeqAnnotator
import pytest

@pytest.mark.skip(reason="test performed manually, not by CI")
def test_fullpipeline(): #{
    # make a copy of the environment
    env = dict(os.environ)

    base_dir = './'

    print('starting the assembly pipeline')


    pre_processor = PreProcessor(base_dir+'src/gsap/', 'AL009126_v11',
                                 base_dir+'tests/data/sample1/SP2_S6_L001_R1_001.fastq.gz',
                                 base_dir+'tests/data/sample1/SP2_S6_L001_R2_001.fastq.gz',
                                 base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.fasta',
                                 base_dir+'tests/data/test_ref_data/',
                                 base_dir+'tests/data/tmp/')


    variant_caller = VariantCaller(base_dir+'src/gsap/', base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.fasta',
                                   base_dir+'tests/data/test_ref_data/addrg_reads.bam',
                                   base_dir+'tests/data/test_ref_data/',
                                   base_dir+'tests/data/tmp/')


    print('######## pre-processing raw reads... ########')
    # run the pre-processing pipeline

    print('quality-based trimming')
    pre_processor.trimSeqReads()

    print('aligning PE reads to reference...')
    pre_processor.alignPairedEndReads('aligned_reads.sam')


    print('sorting aligned reads...')
    pre_processor.sortSAMIntoBAM(base_dir+'tests/data/test_ref_data/aligned_reads.sam', 'sorted_reads.bam') # this command requires lots of memory, 4096MB used

    print('build bam index for sorted_reads.bam')
    pre_processor.buildBamIndex(base_dir+'tests/data/test_ref_data/sorted_reads.bam', base_dir+'tests/data/test_ref_data/sorted_reads.bai')

    print('removing unmapped reads...')
    pre_processor.removeUnmappedReads(base_dir+'tests/data/test_ref_data/sorted_reads.bam', base_dir+'tests/data/test_ref_data/filtered_sample_aligned.bam', base_dir+'tests/data/test_ref_data/unmapped_reads.bam')

    print('build bam index for filtered_sample_aligned.bam')
    pre_processor.buildBamIndex(base_dir+'tests/data/test_ref_data/filtered_sample_aligned.bam', base_dir+'tests/data/test_ref_data/filtered_sample_aligned.bai')

    print('sorting the filtered bam...')
    pre_processor.sortBAMIntoBAM('tests/data/test_ref_data/filtered_sample_aligned.bam', 'filtered_sample_aligned_sorted.bam')
    print('marking duplicate reads...')
    pre_processor.markDuplicates(base_dir+'tests/data/test_ref_data/filtered_sample_aligned_sorted.bam', 'dedup_reads.bam', 'metrics.txt')
    print('adding read group...')
    pre_processor.addReadGroup(base_dir+'tests/data/test_ref_data/dedup_reads.bam', 'addrg_reads.bam')
    print('building bam index from read groups...')
    pre_processor.buildBamIndex(base_dir+'tests/data/test_ref_data/addrg_reads.bam', base_dir+'tests/data/test_ref_data/addrg_reads.bai')
    print('building reference fasta index...')
    pre_processor.buildRefFastaIndex(base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.fasta',
                                     base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.dict')

    print('calling variants...')
    variant_caller.callVariants()

    print('finding a consensus fasta...')
    subprocess.call(['./src/gsap/getconsensusfasta3.sh', base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.fasta', base_dir+'tests/data/test_ref_data/raw_variants.vcf.gz', base_dir+'tests/data/test_ref_data/consensus.fasta', base_dir+'src/gsap'])

    refseq_file = base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.gb'
    refseq_format = 'gb'
    conseq_file = base_dir+'tests/data/test_ref_data/consensus.fasta'
    outdir_base = base_dir+'tests/data/test_ref_data/'
    seqannotator = SeqAnnotator(refseq_file, conseq_file, outdir_base, ref_format=refseq_format)

    in_seq = seqannotator.parseSeqRecordFile(conseq_file)
    ref_seq = seqannotator.parseSeqRecordFile(refseq_file)
    print(f'transfering annotations...')
    in_seq_records = seqannotator.transferAnnotations(in_seq, ref_seq)
    in_seq_records = seqannotator.fillSeqGaps(in_seq_records, ref_seq)
    in_vcf_gz_file = outdir_base+'raw_variants.vcf.gz'
    in_vcf_file = outdir_base+'raw_variants.vcf'
    sh.gunzip('-kf', in_vcf_gz_file)
    in_seq_records = seqannotator.transferVCFAnnotations(in_seq_records, ref_seq, in_vcf_file)
    print(f'transferVCFAnnotations: returned')
    in_seq_records = seqannotator.addHeaderAnnotations(in_seq_records)
    print('writing genbank file...')

    out_gb_file = 'BSB1_annotated.gb'
    seqannotator.writeEMBLFile(in_seq_records, out_gb_file, 'gb')

    print('finished assembling the genome sequence')

    return True
