"""gsap.gsap: provides entry point main()."""

import gzip
import re
import shutil
import subprocess
import uuid
from argparse import ArgumentParser

from gsap import HeaderAnnotationType, PreProcessor, SeqAnnotator, VariantCaller

__version__ = "0.1.0"


def get_parser() -> ArgumentParser:
    """Creates a new argparse argument parser."""
    parser = ArgumentParser("gsap_cli")
    version = "%(prog)s " + __version__
    parser.add_argument("--version", "-v", action="version", version=version)
    parser.add_argument(
        "--fwd_seq",
        "-F",
        help="specify a path to a forward pair-ended short-read sequence file: \
                -F data1.fastq.gz",
        required=True,
    )
    parser.add_argument(
        "--rev_seq",
        "-R",
        help="specify a path to a reverse pair-ended short-read sequence file: \
                -R data2.fastq.gz",
        required=True,
    )
    parser.add_argument(
        "--ref_name",
        "-N",
        help="specify the reference sequence name: -N somename",
        required=True,
    )
    parser.add_argument(
        "--ref_seq_a",
        "-A",
        help="specify a path to a reference genome sequence file in FASTA format: \
                -A data3.fasta",
        required=True,
    )
    parser.add_argument(
        "--ref_seq_b",
        "-B",
        help="specify a path to a reference genome sequence file in Genbank format:\
                -B data3.gb",
        required=True,
    )
    parser.add_argument(
        "--comment",
        "-C",
        help='specify comments for annotating the output Genbank file: \
                -C "comment string"',
        required=True,
    )
    parser.add_argument(
        "--organism",
        "-o",
        help='specify the name of sequenced organism for annotating the output \
                Genbank file: e.g. -o "Bacillus subtilis 168"',
        required=True,
    )
    parser.add_argument(
        "--molecule_type",
        "-m",
        help='specify the sequence molecule type: e.g. -m "DNA"',
        required=True,
    )
    parser.add_argument(
        "--tmp_dir",
        "-T",
        help="specify a path to the directory to hold temp intermediary files: \
                -T ./data_output/",
        required=True,
    )
    parser.add_argument(
        "--out_filepath",
        "-O",
        help="specify a path to the final output file in Genbank format: \
                -O data4.gb",
        required=True,
    )
    return parser


def main(sysargs: list[str] | None = None) -> None:
    """
    The main command line entry point for the GSAP pipeline.

    Args:
    ----
        sysargs : list
            A of arguments as if they were input in the command line.
            when None is given, use sys.argv as an alternative means to parse args.

    """
    sessionid = uuid.uuid4().hex

    parser = get_parser()
    args = parser.parse_args(sysargs)

    base_dir = "./"
    pre_processor = PreProcessor(
        "./src/gsap/",
        args.ref_name,
        args.fwd_seq,
        args.rev_seq,
        args.ref_seq_a,
        args.tmp_dir,
        args.tmp_dir,
    )

    variant_caller = VariantCaller(
        "./src/gsap/",
        args.ref_seq_a,
        f"{args.tmp_dir}/{sessionid}_addrg_reads.bam",
        args.tmp_dir,
        args.tmp_dir,
    )

    print("######## pre-processing raw reads... ########")
    # run the pre-processing pipeline

    print("quality-based trimming")
    pre_processor.trimSeqReads()

    print("aligning PE reads to reference...")
    pre_processor.alignPairedEndReads(f"{sessionid}_aligned_reads.sam")

    print("sorting aligned reads...")
    pre_processor.sortSAMIntoBAM(
        f"{args.tmp_dir}{sessionid}_aligned_reads.sam",
        f"{sessionid}_sorted_reads.bam",
    )  # this command requires lots of memory, 4096MB used

    print("build bam index for sorted_reads.bam")
    pre_processor.buildBamIndex(
        f"{args.tmp_dir}{sessionid}_sorted_reads.bam",
        f"{args.tmp_dir}{sessionid}_sorted_reads.bai",
    )

    print("removing unmapped reads...")
    pre_processor.removeUnmappedReads(
        f"{args.tmp_dir}{sessionid}_sorted_reads.bam",
        f"{args.tmp_dir}{sessionid}_filtered_sample_aligned.bam",
        f"{args.tmp_dir}{sessionid}_unmapped_reads.bam",
    )

    print("build bam index for filtered_sample_aligned.bam")
    pre_processor.buildBamIndex(
        f"{args.tmp_dir}{sessionid}_filtered_sample_aligned.bam",
        f"{args.tmp_dir}{sessionid}_filtered_sample_aligned.bai",
    )

    print("sorting the filtered bam...")
    pre_processor.sortBAMIntoBAM(
        f"{args.tmp_dir}{sessionid}_filtered_sample_aligned.bam",
        f"{sessionid}_filtered_sample_aligned_sorted.bam",
    )
    print("marking duplicate reads...")
    pre_processor.markDuplicates(
        f"{args.tmp_dir}{sessionid}_filtered_sample_aligned_sorted.bam",
        f"{sessionid}_dedup_reads.bam",
        f"{sessionid}_metrics.txt",
    )
    print("adding read group...")
    pre_processor.addReadGroup(
        f"{args.tmp_dir}{sessionid}_dedup_reads.bam",
        f"{sessionid}_addrg_reads.bam",
    )
    print("building bam index from read groups...")
    pre_processor.buildBamIndex(
        f"{args.tmp_dir}{sessionid}_addrg_reads.bam",
        f"{args.tmp_dir}{sessionid}_addrg_reads.bai",
    )
    print("building reference fasta index...")
    re_match = re.search("(.*)\\.(gb|embl|fasta)$", args.ref_seq_a)
    if type(re_match) is not re.Match:
        raise Exception("The ref_seq_a file has unsupported extension")
    seq_filepath_noext = re_match.groups()[0]
    pre_processor.buildRefFastaIndex(args.ref_seq_a, f"{seq_filepath_noext}.dict")

    print("calling variants...")
    variant_caller.callVariants(f"{sessionid}_raw_variants.vcf.gz")

    print("finding a consensus fasta...")
    subprocess.call(
        [
            "./src/gsap/getconsensusfasta3.sh",
            args.ref_seq_a,
            f"{args.tmp_dir}{sessionid}_raw_variants.vcf.gz",
            f"{args.tmp_dir}{sessionid}_consensus.fasta",
            base_dir + "src/gsap",
        ]
    )

    refseq_format = "gb"
    conseq_file = f"{args.tmp_dir}{sessionid}_consensus.fasta"
    outdir_base = args.tmp_dir
    seqannotator = SeqAnnotator(
        args.ref_seq_b, conseq_file, outdir_base, ref_format=refseq_format
    )

    in_seq = seqannotator.parseSeqRecordFile(conseq_file)
    ref_seq = seqannotator.parseSeqRecordFile(args.ref_seq_b)
    print("transferring annotations...")
    in_seq_records = seqannotator.transferAnnotations(in_seq, ref_seq)
    in_seq_records = seqannotator.fillSeqGaps(in_seq_records, ref_seq)
    in_vcf_gz_file = outdir_base + f"{sessionid}_raw_variants.vcf.gz"
    in_vcf_file = outdir_base + f"{sessionid}_raw_variants.vcf"

    with gzip.open(in_vcf_gz_file, "rb") as f_in:
        with open(in_vcf_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    in_seq_records = seqannotator.transferVCFAnnotations(
        in_seq_records, ref_seq, in_vcf_file
    )
    print("transferVCFAnnotations: returned")
    # add header annotations

    header_annotations: HeaderAnnotationType = {}
    header_annotations["comment"] = str(args.comment)
    header_annotations["organism"] = str(args.organism)
    header_annotations["taxonomy"] = "taxon:unknown"  # args.taxonomy

    header_annotations["molecule_type"] = str(args.molecule_type)
    in_seq_records = seqannotator.addHeaderAnnotations(
        in_seq_records, header_annotations
    )
    print("writing genbank file...")

    out_gb_file = f"{sessionid}_output.gb"
    seqannotator.writeEMBLFile(in_seq_records, out_gb_file, "gb")

    shutil.copyfile(outdir_base + out_gb_file, args.out_filepath)

    print("finished assembling the genome sequence")
