import os
from pytest import raises

from preprocessing import PreProcessor

base_dir = './data/'
testdata_dir = './test_data/'
output_dir = './test_data_wd/'

pre_processor = PreProcessor(base_dir, 'AL009126_v10',
                             base_dir+'sample1/SP2_S6_L001_R1_001.fastq.gz',
                             base_dir+'sample1/SP2_S6_L001_R2_001.fastq.gz',
                             base_dir+'RefGenome/AL009126_v10/AL009126.fasta',
                             output_dir)


def test_alignPairedEndReads():
    pre_processor.alignPairedEndReads('aligned_reads.sam')
    assert os.path.exists(output_dir+'aligned_reads.sam') == True

def test_sortSAMIntoBAM():
    pre_processor.sortSAMIntoBAM(output_dir+'aligned_reads.sam', 'sorted_reads.bam')
    assert os.path.exists(output_dir+'sorted_reads.bam') == True

def test_removeUnmappedReads():
    pre_processor.removeUnmappedReads(output_dir+'sorted_reads.bam', output_dir+'filtered_sample_aligned.bam', output_dir+'unmapped_reads.bam')

    assert os.path.exists(output_dir+'filtered_sample_aligned.bam') == True
    assert os.path.exists(output_dir+'unmapped_reads.bam') == True

def test_buildBamIndex():
    pre_processor.buildBamIndex(output_dir+'sorted_reads.bam')
    assert os.path.exists(output_dir+'sorted_reads.bai') == True

def test_buildRefFastaIndex():
    pre_processor.buildRefFastaIndex('./data/RefGenome/AL009126_EVOt2/AL009126_evot2.fasta', './data/RefGenome/AL009126_EVOt2/AL009126_evot2.dict')
