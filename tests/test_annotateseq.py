__author__ = "spark"

import gzip
import os
import shutil
import subprocess

import pytest

from gsap import SeqAnnotator

base_dir = "./src/gsap/"
testdata_dir = "./tests/data/test_ref_data/"
output_dir = "test_data_wd"

# arrange fixtures
@pytest.fixture
def gsap_test_common(tmp_path):
    src_base_dir = "./src/gsap/"
    test_ref_data_dir = "./tests/data/test_ref_data/"
    tmp_test_data_dir = tmp_path / "test_data_wd"
    tmp_test_data_dir.mkdir()

    refseq_file = "./tests/data/RefGenome/AL009126_v11/AL009126.gb"
    refseq_format = "gb"
    conseq_file = test_ref_data_dir+"consensus.fasta"
    seqannotator = SeqAnnotator(
        refseq_file, conseq_file, str(tmp_test_data_dir), ref_format=refseq_format
    )

    def pytest_runtest_teardown():
        tmp_test_data_dir.remove()

    return (
        seqannotator,
        test_ref_data_dir,
        str(tmp_test_data_dir) + "/",
        pytest_runtest_teardown,
    )

def test_getconsensusfasta3(gsap_test_common):
    seqannotator, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    subprocess.call(
        [
            "./src/gsap/getconsensusfasta3.sh",
            "./tests/data/RefGenome/AL009126_v11/AL009126.fasta",
            test_ref_data_dir+"raw_variants.vcf.gz",
            tmp_test_data_dir+"consensus.fasta",
            base_dir + "src/gsap",
        ],
        shell=False
    )

    assert os.path.exists(tmp_test_data_dir + "consensus.fasta") == True


def test_annotateseq(gsap_test_common):
    seqannotator, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common

    in_seq = seqannotator.parseSeqRecordFile(seqannotator._conseq_file)
    ref_seq = seqannotator.parseSeqRecordFile(seqannotator._refseq_file)
    print("transferring annotations...")
    in_seq_records = seqannotator.transferAnnotations(in_seq, ref_seq)
    in_seq_records = seqannotator.fillSeqGaps(in_seq_records, ref_seq)
    in_vcf_gz_file = test_ref_data_dir + "raw_variants.vcf.gz"
    in_vcf_file = tmp_test_data_dir + "raw_variants.vcf"

    with gzip.open(in_vcf_gz_file, 'rb') as f_in:
      with open(in_vcf_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    in_seq_records = seqannotator.transferVCFAnnotations(
        in_seq_records, ref_seq, in_vcf_file
    )
    print("transferVCFAnnotations: returned")
    header_annotations = {}
    header_annotations["comment"] = (
        "CBCB's BSB1 strain of B. subtilis: sequenced using Illumina and assembled against the reference strain 168 (AL009126.3) dated 26-Feb-2014. The assembly was done using bowtie-2, and the variant analysis was done using GATK (Genome Analysis ToolKit). The annotations (gene, CDS, tRNA and rRNA) from the reference strain AL009126 were transferred using RATT (Rapid Annotation Transfer Tool). The VCF file from the GATK pipeline was used to add annotations on variant calls."
    )
    header_annotations["organism"] = "Bacillus subtilis subsp. BSB1 from CBCB"
    header_annotations["taxonomy"] = [
        "Bacteria",
        "Firmicutes",
        "Bacilli",
        "Bacillales",
        "Bacillaceae",
        "Bacillus",
    ]
    header_annotations["molecule_type"] = "DNA"
    in_seq_records = seqannotator.addHeaderAnnotations(
        in_seq_records, header_annotations
    )
    print("writing genbank file...")

    out_gb_file = "BSB1_annotated.gb"
    seqannotator.writeEMBLFile(in_seq_records, out_gb_file, "gb")

    assert os.path.exists(tmp_test_data_dir + "BSB1_annotated.gb") == True

def teardown(gsap_test_common):
    seqannotator, test_ref_data_dir, tmp_test_data_dir, teardown = gsap_test_common
    # tear down tmp file and dir
    # os.remove(tmp_input_file_path)
    # os.remove(tmp_output_file_path)
    # os.rmdir(tmp_output_dir)
    teardown()
