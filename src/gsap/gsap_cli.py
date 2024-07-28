"""gsap.gsap: provides entry point main()."""

import argparse
from pathlib import Path

import shutil

import uuid

from gsap import PreProcessor, SeqAnnotator, VariantCaller

__version__ = "0.1.0"

def get_parser():
    """Creates a new argparse argument parser."""
    parser = argparse.ArgumentParser('gsap_cli')
    version = '%(prog)s ' + __version__
    parser.add_argument('--version', '-v', action='version', version=version)
    parser.add_argument('--fwd_seq', '-F', 
                        help='specify a path to a forward pair-ended short-read sequence file: -F data1.fastq.gz', required=True )
    parser.add_argument('--rev_seq', '-R', 
                        help='specify a path to a reverse pair-ended short-read sequence file: -R data2.fastq.gz', required=True)
    parser.add_argument('--ref_name', '-N', 
                        help='specify the reference sequence name: -N somename', required=True)
    parser.add_argument('--ref_seq_a', '-A', 
                        help='specify a path to a reference genome sequence file in FASTA format: -A data3.fasta', required=True)
    parser.add_argument('--ref_seq_b', '-B', 
                        help='specify a path to a reference genome sequence file in Genbank format: -B data3.gb', required=True)
    parser.add_argument('--comment', '-C', 
                        help='specify comments for annotating the output Genbank file: -C "comment string"', required=True)
    parser.add_argument('--organism', '-o', 
                        help='specify the name of sequenced organism for annotating the output Genbank file: e.g. -o "Bacillus subtilis 168"', required=True)
    parser.add_argument('--molecule_type', '-m', 
                        help='specify the sequence molecule type: e.g. -m "DNA"', required=True)
    parser.add_argument('--out_filepath', '-O', 
                        help='specify a path to the final output file in Genbank format: -O data4.gb', required=True)
    return parser


def main(args : list[str] | None=None):
    """
    The main command line entry point for the GSAP pipeline.
    Args:
        args : list
            A of arguments as if they were input in the command line. 
            when None is given, use sys.argv as an alternative means to parse args.
    """

    sessionid = uuid.uuid4().hex

    parser = get_parser()
    args = parser.parse_args(args)

    base_dir = "./"
    pre_processor = PreProcessor("./src/gsap/",
                                 args.ref_name,
                                 args.fwd_seq,
                                 args.rev_seq,
                                 args.ref_seq_a,
                                 "./tests/data/test_ref_data/")

    variant_caller = VariantCaller("./src/gsap/",
                                   args.ref_seq_a,
                                   f"{base_dir}/tests/data/test_ref_data/{sessionid}_addrg_reads.bam",
                                   f"{base_dir}/tests/test_ref_data/")

    print("######## pre-processing raw reads... ########")
    # run the pre-processing pipeline

    print("quality-based trimming")
    pre_processor.trimSeqReads()

    print("aligning PE reads to reference...")
    pre_processor.alignPairedEndReads(f"{sessionid}_aligned_reads.sam")

    print("sorting aligned reads...")
    pre_processor.sortSAMIntoBAM(
        base_dir + f"tests/data/test_ref_data/{sessionid}_aligned_reads.sam", f"{sessionid}_sorted_reads.bam"
    )  # this command requires lots of memory, 4096MB used

    print("build bam index for sorted_reads.bam")
    pre_processor.buildBamIndex(
        base_dir + f"tests/data/test_ref_data/{sessionid}_sorted_reads.bam",
        base_dir + f"tests/data/test_ref_data/{sessionid}_sorted_reads.bai",
    )

    print("removing unmapped reads...")
    pre_processor.removeUnmappedReads(
        base_dir + f"tests/data/test_ref_data/{sessionid}_sorted_reads.bam",
        base_dir + f"tests/data/test_ref_data/{sessionid}_filtered_sample_aligned.bam",
        base_dir + f"tests/data/test_ref_data/{sessionid}_unmapped_reads.bam",
    )

    print("build bam index for filtered_sample_aligned.bam")
    pre_processor.buildBamIndex(
        base_dir + f"tests/data/test_ref_data/{sessionid}_filtered_sample_aligned.bam",
        base_dir + f"tests/data/test_ref_data/{sessionid}_filtered_sample_aligned.bai",
    )

    print("sorting the filtered bam...")
    pre_processor.sortBAMIntoBAM(
        f"tests/data/test_ref_data/{sessionid}_filtered_sample_aligned.bam",
        f"{sessionid}_filtered_sample_aligned_sorted.bam",
    )
    print("marking duplicate reads...")
    pre_processor.markDuplicates(
        base_dir + f"tests/data/test_ref_data/{sessionid}_filtered_sample_aligned_sorted.bam",
        f"{sessionid}_dedup_reads.bam",
        f"{sessionid}_metrics.txt",
    )
    print("adding read group...")
    pre_processor.addReadGroup(
        base_dir + f"tests/data/test_ref_data/{sessionid}_dedup_reads.bam", f"{sessionid}_addrg_reads.bam"
    )
    print("building bam index from read groups...")
    pre_processor.buildBamIndex(
        base_dir + f"tests/data/test_ref_data/{sessionid}_addrg_reads.bam",
        base_dir + f"tests/data/test_ref_data/{sessionid}_addrg_reads.bai",
    )
    print("building reference fasta index...")
    pre_processor.buildRefFastaIndex(
        args.ref_seq_a,
        f"{args.ref_seq_a}.{sessionid}.dict"
    )

    print("calling variants...")
    variant_caller.callVariants()

    print("finding a consensus fasta...")
    subprocess.call(
        [
            "./src/gsap/getconsensusfasta3.sh",
            args.ref_seq_a,
            base_dir + "tests/data/test_ref_data/raw_variants.vcf.gz",
            base_dir + f"tests/data/test_ref_data/{sessionid}_consensus.fasta",
            base_dir + "src/gsap",
        ]
    )

    refseq_format = "gb"
    conseq_file = base_dir + f"tests/data/test_ref_data/{sessionid}_consensus.fasta"
    outdir_base = base_dir + "tests/data/test_ref_data/"
    seqannotator = SeqAnnotator(
        args.ref_seq_b, conseq_file, outdir_base, ref_format=refseq_format
    )

    in_seq = seqannotator.parseSeqRecordFile(conseq_file)
    ref_seq = seqannotator.parseSeqRecordFile(refseq_file)
    print("transferring annotations...")
    in_seq_records = seqannotator.transferAnnotations(in_seq, ref_seq)
    in_seq_records = seqannotator.fillSeqGaps(in_seq_records, ref_seq)
    in_vcf_gz_file = outdir_base + "raw_variants.vcf.gz"
    in_vcf_file = outdir_base + f"{sessionid}_raw_variants.vcf"
    sh.gunzip("-kf", in_vcf_gz_file)
    in_seq_records = seqannotator.transferVCFAnnotations(
        in_seq_records, ref_seq, in_vcf_file
    )
    print("transferVCFAnnotations: returned")
    # add header annotations

    header_annotations = {}
    header_annotations["comment"] = args.comment
    header_annotations["organism"] = args.organism
    #header_annotations["taxonomy"] = args.taxonomy

    header_annotations["molecule_type"] = args.molecule_type
    in_seq_records = seqannotator.addHeaderAnnotations(
        in_seq_records, header_annotations
    )
    print("writing genbank file...")

    out_gb_file = f"{sessionid}_output.gb"
    seqannotator.writeEMBLFile(in_seq_records, out_gb_file, "gb")


    shutil.copyfile(
        out_gb_file,
        args.out_filepath
    )

    print("finished assembling the genome sequence")

