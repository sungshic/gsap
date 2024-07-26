__author__ = 'spark'

import os
import shutil
import pytest
from pytest import raises
from gsap import VariantCaller, PreProcessor

base_dir = './src/gsap/'
testdata_dir = './tests/data/test_ref_data/'
output_dir = 'test_data_wd'

# arrange fixtures
@pytest.fixture
def gsap_test_common(tmp_path):
    src_base_dir = './src/gsap/'
    test_ref_data_dir = './tests/data/test_ref_data/'
    tmp_test_data_dir = tmp_path / 'test_data_wd'
    tmp_test_data_dir.mkdir()

    shutil.copyfile(test_ref_data_dir+'addrg_reads.bam', str(tmp_test_data_dir)+'/addrg_reads.bam')
    variant_caller = VariantCaller(src_base_dir, './tests/data/RefGenome/AL009126_v11/AL009126.fasta',
                                   str(tmp_test_data_dir) + '/addrg_reads.bam',
                                   str(tmp_test_data_dir)+'/',
                                   str(tmp_test_data_dir)+'/')
    def pytest_runtest_teardown():
        tmp_test_data_dir.remove()

    return variant_caller, test_ref_data_dir, str(tmp_test_data_dir)+"/", pytest_runtest_teardown

def test_callVariants(gsap_test_common):
    variant_caller, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    variant_caller.callVariants("raw_variants.vcf.gz")
    assert os.path.exists(tmp_test_data_dir+"raw_variants.vcf.gz") == True


def teardown(gsap_test_common):
    variant_caller, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    # tear down tmp file and dir
    #os.remove(tmp_input_file_path)
    #os.remove(tmp_output_file_path)
    #os.rmdir(tmp_output_dir)
    teardown()
