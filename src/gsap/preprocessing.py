"""
Implements functions for analysing and assembling genome DNA data.

Functions perform pre-processing steps to transform raw short-read sequence data,
initially in FASTQ format into other data forms as required by tools such as
samtools, picard-tools, and bowtie2.

Author: Sunny Park

"""

import subprocess

__all__ = [
    "PreProcessor",
]


class PreProcessor:
    """
    Offers steps for assembling a genome sequence from raw short read sequences.

    This class is used to manage a variety of toolset functions employed by GSAP
    to accomplish various sequence assembly-related tasks. The tasks defined as
    class methods utilize external functions implemented in Python and Java,
    as well as other bioinformatic tools distributed as Linux command line binaries.

    Following are the details of constructor parameters to initialise
    key input data for GSAP's sequence assembly operations

    :param home_dir: the directory location of the source files
    :type home_dir: str
    :param refgenome_name: name of a reference genome
    :type refgenome_name: str
    :param rawdata_filepath: the file path of a gzipped FASTQ file containing
        forward strands of raw paired-end short read sequences from NGS devices
        like Illumina. In case of rawdata_paired2_filepath being None, this
        parameter is for a single-read sequence file in gzipped FASTQ format.
    :type rawdata_filepath: str
    :param rawdata_paired2_filepath: the file path of a gzipped FASTQ file
        containing reverse strands of raw paired-end short read sequences
        from NGS devices like Illumina
    :type rawdata_paired2_filepath: str | None
    :param refgenome_filepath: the file path of the reference genome in FASTA
        format.
    :type refgenome_filepath: str
    :param data_outdir: the directory where intermediary files, which may
        or may not persist after the end of assembly, are stored during the
        assembly operation.
    :type data_outdir: str, optional
    :param tmpdir: the directory where intermediary transient files, if any,
        are stored during the assembly operation. [Note: currently not in use]
    :type tmpdir: str, optional
    """

    def __init__(
        self,
        home_dir: str,
        refgenome_name: str,
        rawdata_filepath: str,
        rawdata_paired2_filepath: str | None = None,
        refgenome_filepath: str | None = None,
        data_outdir: str = "data/out/",
        tmpdir: str = "data/tmp",
    ):
        """Constructor method."""
        self._rawdata_filepath = rawdata_filepath
        self._home_dir = home_dir
        if rawdata_paired2_filepath is not None:
            self._isPairedEnd = True
            self._rawdata_paired2_filepath = rawdata_paired2_filepath
        if refgenome_filepath is not None:
            self._refgenome_filepath = refgenome_filepath
            self.indexRefGenome(self._refgenome_filepath, refgenome_name)
        self._refgenome_name = refgenome_name
        self._data_outdir = data_outdir
        self._tmpdir = tmpdir

    def qualityCheck(self, fastq_filepath: str) -> None:
        """
        Quality Check method using the command line tool fastqc.

        :param fastq_filepath: the path to locate a FASTQ file
        :type fastq_filepath: str
        """
        subprocess.call(["/usr/bin/fastqc", fastq_filepath], shell=False)

    def trimSeqReads(
        self,
        leading: int = 28,
        trailing: int = 28,
        headcrop: int = 15,
        minlen: int = 36,
    ) -> None:
        """
        A method to trim sequences based on Quality (Phred) scores, with base 33.

        :param leading: Cut bases off the start of a read, if below
            a threshold quality, default used: 28
        :type leading: int
        :param trailing: Cut bases off the end of a read, if below
            a threshold quality, default used: 28
        :type trailing: int
        :param headcrop: specified number of bases cur from the start of
            the read, default used: 15
        :type headcrop: int
        :param minlen: drop the read if it is below a specified length,
            default used: 36
        :type minlen: int
        """
        fastq_filepath_for = self._rawdata_filepath
        fastq_filepath_rev = self._rawdata_paired2_filepath
        out_forward_paired = self._data_outdir + "out_for_paired.fq.gz"
        out_forward_unpaired = self._data_outdir + "out_for_unpaired.fq.gz"
        out_rev_paired = self._data_outdir + "out_rev_paired.fq.gz"
        out_rev_unpaired = self._data_outdir + "out_rev_unpaired.fq.gz"

        try:
            subprocess.call(
                [
                    "/usr/bin/java",
                    "-jar",
                    "/usr/share/java/trimmomatic.jar",
                    "PE",
                    "-phred33",
                    fastq_filepath_for,
                    fastq_filepath_rev,
                    out_forward_paired,
                    out_forward_unpaired,
                    out_rev_paired,
                    out_rev_unpaired,
                    "LEADING:" + str(leading),
                    "TRAILING:" + str(trailing),
                    "HEADCROP:" + str(headcrop),
                    "MINLEN:" + str(minlen),
                ],
                shell=False,
            )  # , sort_order])
            # replace the raw files with trimmed sequences
            self._rawdata_filepath = out_forward_paired
            self._rawdata_paired2_filepath = out_rev_paired
        except Exception as e:
            print(e)

    def assembleDeNovo(self) -> None:  # , kmer_range):
        """Assembles raw short-read sequences de novo without a ref genome."""
        out_forward_paired = self._data_outdir + "out_for_paired.fq.gz"
        out_rev_paired = self._data_outdir + "out_rev_paired.fq.gz"
        subprocess.call(
            [
                "/gsap/src/gsap/toolset/SPAdes-4.0.0-Linux/bin/spades.py",
                "--pe1-1",
                out_forward_paired,
                "--pe1-2",
                out_rev_paired,
                "-o",
                self._data_outdir + "denovo/spades/output",
            ],
            shell=False,
        )

    def indexRefGenome(self, refgenome_filepath: str, refgenome_name: str) -> None:
        """Builds a Bowtie index from a reference DNA sequence."""
        subprocess.call(
            ["/usr/bin/bowtie2-build", refgenome_filepath, refgenome_name], shell=False
        )

    def alignPairedEndReads(self, output_filename: str) -> None:
        """Aligns paired-end short-read sequences w.r.t. a ref genome."""
        self._alignPairedEndReads(
            output_filename,
            self._refgenome_name,
            self._rawdata_filepath,
            self._rawdata_paired2_filepath,
        )

    def _alignPairedEndReads(
        self,
        output_filename: str,
        ref_name: str,
        reads1_filepath: str,
        reads2_filepath: str,
    ) -> None:
        """Private method to invoke bowtie2 for alignment w.r.t. a ref genome."""
        subprocess.call(
            [
                "/usr/bin/bowtie2",
                "-x",
                ref_name,
                "-1",
                reads1_filepath,
                "-2",
                reads2_filepath,
                "-S",
                self._data_outdir + output_filename,
            ],
            shell=False,
        )

    def sortSAMIntoBAM(
        self, insam_filepath: str, outbam_filename: str, sort_order: str = "coordinate"
    ) -> None:
        """Sorts a SAM file and outputs the sorted sequence as a BAM file."""
        try:
            subprocess.call(
                [
                    "/usr/bin/picard-tools",
                    "SortSam",
                    "INPUT=",
                    insam_filepath,
                    "OUTPUT=",
                    self._data_outdir + outbam_filename,
                    "SORT_ORDER=coordinate",
                ],
                shell=False,
            )  # , sort_order])
        except Exception as e:
            print(e)

    def sortBAMIntoBAM(
        self, inbam_filepath: str, outbam_filename: str, sort_order: str = "coordinate"
    ) -> None:
        """Sorts a BAM file and outputs the sorted sequence as a BAM file."""
        try:
            subprocess.call(
                [
                    "/usr/bin/picard-tools",
                    "SortSam",
                    "INPUT=",
                    inbam_filepath,
                    "OUTPUT=",
                    self._data_outdir + outbam_filename,
                    "SORT_ORDER=coordinate",
                ],
                shell=False,
            )  # , sort_order])
        except Exception as e:
            print(e)

    def markDuplicates(
        self, inbam_filepath: str, outbam_filename: str, outmetrics_filename: str
    ) -> None:
        """Locates and tags duplicate reads in a BAM or SAM file."""
        try:
            subprocess.call(
                [
                    "/usr/bin/picard-tools",
                    "MarkDuplicates",
                    "INPUT=",
                    inbam_filepath,
                    "OUTPUT=",
                    self._data_outdir + outbam_filename,
                    "METRICS_FILE=",
                    self._data_outdir + outmetrics_filename,
                ],
                shell=False,
            )
        except Exception as e:
            print(e)

    def addReadGroup(self, inbam_filepath: str, outbam_filename: str) -> None:
        """Assigns all the reads in a file to a single new read-group."""
        try:
            subprocess.call(
                [
                    "/usr/bin/picard-tools",
                    "AddOrReplaceReadGroups",
                    "INPUT=",
                    inbam_filepath,
                    "OUTPUT=",
                    self._data_outdir + outbam_filename,
                    "RGID=",
                    "group1",
                    "RGLB=",
                    "lib1",
                    "RGPL=",
                    "illumina",
                    "RGPU=",
                    "unit1",
                    "RGSM=",
                    "sample1",
                ],
                shell=False,
            )
        except Exception as e:
            print(e)

    def buildBamIndex(self, inbam_filepath: str, outbai_filepath: str) -> None:
        """Generate a BAM index (.bai) file from a BAM (.bam) file."""
        try:
            subprocess.call(
                [
                    "/usr/bin/picard-tools",
                    "BuildBamIndex",
                    "--INPUT",
                    inbam_filepath,
                    "--OUTPUT",
                    outbai_filepath,
                ],
                shell=False,
            )
        except Exception as e:
            print(e)

    def buildRefFastaIndex(
        self, inreffasta_filepath: str, outdict_filepath: str
    ) -> None:
        """Create a SAM/BAM file and a faidx file from a reference sequence in FASTA."""
        try:
            subprocess.call(["/usr/bin/rm", outdict_filepath], shell=False)
            # create a dictionary file from fasta
            subprocess.call(
                [
                    "/usr/bin/picard-tools",
                    "CreateSequenceDictionary",
                    "R=",
                    inreffasta_filepath,
                    "O=",
                    outdict_filepath,
                ],
                shell=False,
            )
            # create a fasta index file
            subprocess.call(
                ["/usr/bin/samtools", "faidx", inreffasta_filepath], shell=False
            )
        except Exception as e:
            print(e)

    def removeUnmappedReads(
        self,
        inbam_filepath: str,
        mapped_outbam_filepath: str,
        unmapped_outbam_filepath: str,
    ) -> None:
        """Splits sequence segments in BAM into mapped and unmapped sequence files."""
        try:
            with open(mapped_outbam_filepath, "w") as f_out:
                subprocess.call(
                    ["/usr/bin/samtools", "view", "-hbF", "4", inbam_filepath],
                    shell=False,
                    stdout=f_out,
                )
            with open(unmapped_outbam_filepath, "w") as f_out:
                subprocess.call(
                    ["/usr/bin/samtools", "view", "-hbf", "4", inbam_filepath],
                    shell=False,
                    stdout=f_out,
                )
        except Exception as e:
            print(e)
