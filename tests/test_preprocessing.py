__author__ = 'spark'

import os
import shutil
import pytest
from pytest import raises
from gsap import VariantCaller, PreProcessor

# arrange fixtures
@pytest.fixture
def gsap_test_common(tmp_path):
    src_base_dir = './src/gsap/'
    test_ref_data_dir = './tests/data/test_ref_data/'
    tmp_test_data_dir = tmp_path / 'test_data_wd'
    tmp_test_data_dir.mkdir()
    pre_processor = PreProcessor(src_base_dir, 'AL009126_v11',
                                 './tests/data/sample1/SP2_S6_L001_R1_001.fastq.gz',
                                 './tests/data/sample1/SP2_S6_L001_R2_001.fastq.gz',
                                 './tests/data/RefGenome/AL009126_v10/AL009126.fasta',
                                 str(tmp_test_data_dir)+"/")

    def pytest_runtest_teardown():
        tmp_test_data_dir.remove()

    return pre_processor, test_ref_data_dir, str(tmp_test_data_dir)+"/", pytest_runtest_teardown

def test_alignPairedEndReads(gsap_test_common):
    pre_processor, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    print(f"test_alignPairedEndReads: {test_ref_data_dir}, {tmp_test_data_dir}")
    shutil.copyfile(test_ref_data_dir+'addrg_reads.bam', tmp_test_data_dir+'addrg_reads.bam')
    tmp_output_file_path = tmp_test_data_dir + "aligned_reads.sam"
    pre_processor.alignPairedEndReads('aligned_reads.sam')
    assert os.path.exists(tmp_output_file_path) == True

def test_sortSAMIntoBAM(gsap_test_common):
    pre_processor, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    shutil.copyfile(test_ref_data_dir+'aligned_reads.sam', tmp_test_data_dir+'aligned_reads.sam')
    tmp_output_file_path = tmp_test_data_dir + "sorted_reads.bam"
    pre_processor.sortSAMIntoBAM(tmp_test_data_dir+'aligned_reads.sam', 'sorted_reads.bam')
    assert os.path.exists(tmp_output_file_path) == True

def test_removeUnmappedReads(gsap_test_common):
    pre_processor, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    shutil.copyfile(test_ref_data_dir+'sorted_reads.bam', tmp_test_data_dir+'sorted_reads.bam')
    pre_processor.removeUnmappedReads(tmp_test_data_dir+'sorted_reads.bam', tmp_test_data_dir+'filtered_sample_aligned.bam', tmp_test_data_dir+'unmapped_reads.bam')

    assert os.path.exists(tmp_test_data_dir+'filtered_sample_aligned.bam') == True
    assert os.path.exists(tmp_test_data_dir+'unmapped_reads.bam') == True

def test_buildBamIndex(gsap_test_common):
    pre_processor, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    shutil.copyfile(test_ref_data_dir+'sorted_reads.bam', tmp_test_data_dir+'sorted_reads.bam')
    pre_processor.buildBamIndex(tmp_test_data_dir+'sorted_reads.bam', tmp_test_data_dir+'sorted_reads.bai')
    assert os.path.exists(tmp_test_data_dir+'sorted_reads.bai') == True

def test_buildRefFastaIndex(gsap_test_common):
    pre_processor, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    shutil.copyfile(test_ref_data_dir+'sorted_reads.bam', tmp_test_data_dir+'sorted_reads.bam')
    pre_processor.buildRefFastaIndex('./tests/data/RefGenome/AL009126_v11/AL009126.fasta', tmp_test_data_dir+'AL009126.dict')
    assert os.path.exists(tmp_test_data_dir+'AL009126.dict') == True

def teardown(gsap_test_common):
    pre_processor, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    teardown()

