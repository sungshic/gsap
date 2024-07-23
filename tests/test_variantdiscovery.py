__author__ = 'spark'

import os
import shutil
from gsap import VariantCaller

base_dir = './src/gsap/'
testdata_dir = './tests/data/test_ref_data/'
output_dir = 'test_data_wd'

def test_callVariants(tmp_path):
    tmp_input_file_path = tmp_path / "addrg_reads.bam"
    tmp_output_dir = tmp_path / output_dir
    tmp_output_dir.mkdir()
    tmp_output_file_path = tmp_output_dir / "raw_variants.vcf"
    shutil.copyfile(testdata_dir+'addrg_reads.bam', tmp_input_file_path)

    variant_caller = VariantCaller(base_dir, './tests/data/RefGenome/AL009126_v10/AL009126.fasta',
                                   str(tmp_input_file_path),
                                   str(tmp_output_dir)+'/')

    variant_caller.callVariants()
    print(f"tmp_output_file_path: {tmp_output_file_path}")
    assert os.path.exists(tmp_output_file_path) == True

    # tear down tmp file and dir
    os.remove(tmp_input_file_path)
    os.remove(tmp_output_file_path)
    os.rmdir(tmp_output_dir)
