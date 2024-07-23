# copied from http://www.biostars.org/p/6970/
# You can do this with a combination of Biopython for writing the Fasta files and pysam for reading the BAM files:
# Run it with:
# % python bam2fasta.py your_file.bam
# and the output will be in your_file.fa


import os
import sys

import pysam
from Bio import SeqIO, Seq, SeqRecord

def main(in_file):
    out_file = "%s.fa" % os.path.splitext(in_file)[0]
    with open(out_file, "w") as out_handle:
        # Write records from the BAM file one at a time to the output file.
        # Works lazily as BAM sequences are read so will handle large files.
        SeqIO.write(bam_to_rec(in_file), out_handle, "fasta")

def bam_to_rec(in_file):
    """Generator to convert BAM files into Biopython SeqRecords.
    """
    bam_file = pysam.Samfile(in_file, "rb")
    for (read in bam_file):
        seq = Seq.Seq(read.seq)
        if (read.is_reverse):
            seq = seq.reverse_complement()
        rec = SeqRecord.SeqRecord(seq, read.qname, "", "")
        yield rec

if (__name__ == "__main__"):
    main(*sys.argv[1:])
