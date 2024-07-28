"""
Implements glue functions for variant analysis in genome DNA data.

Functions perform variant calling against a reference genome via GATK.

Author: Sunny Park

"""

import subprocess

__all__ = [
    "VariantCaller",
]


class VariantCaller:
    """
    Offers a variant calling method from the sequency assembly input in BAM format.

    Following are the details of constructor parameters to initialise
    key input data for GSAP's variant calling operations

    :param home_dir: the directory location of the source files
    :type home_dir: str
    :param refgenome_filepath: the file path of a reference genome in FASTA
        format.
    :type refgenome_filepath: str
    :param inputbam_filepath: the file path of an assembled genome sequence in
        BAM format.
    :type inputbam_filepath: str
    :param data_outdir: the output directory where a gzipped VCF file carrying SNPs
        and INDELs is stored.
    :type data_outdir: str, optional
    :param tmpdir: the directory where intermediary transient files, if any,
        are stored during the variant calling operation. [Note: currently not in use]
    :type tmpdir: str, optional
    """

    def __init__(
        self,
        home_dir: str,
        refgenome_filepath: str,
        inputbam_filepath: str,
        data_outdir: str = "data/out/",
        tmpdir: str = "data/tmp",
    ):
        """Constructor method."""
        self._refgenome_filepath = refgenome_filepath
        self._inputbam_filepath = inputbam_filepath
        self._data_outdir = data_outdir
        self._tmpdir = tmpdir
        self._home_dir = home_dir

    def callVariants(self, outputvcf_filename: str = "raw_variants.vcf.gz") -> None:
        """A wrapper function to invoke the variant calling operation."""
        try:
            self._callVariants(
                self._refgenome_filepath, self._inputbam_filepath, outputvcf_filename
            )
        except Exception as e:
            print(e)

    def _callVariants(
        self, refgenome_filepath: str, inputbam_filepath: str, outputvcf_filename: str
    ) -> None:
        """Performs variant calling via GATK."""
        try:
            subprocess.call(
                [
                    "/usr/bin/java",
                    "-Xmx4g",
                    "-Djava.io.tmpdir=" + self._tmpdir,
                    "-jar",
                    self._home_dir
                    + "toolset/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar",
                    "HaplotypeCaller",
                    "--ploidy",
                    "1",
                    "-R",
                    refgenome_filepath,
                    "-I",
                    inputbam_filepath,
                    "-O",
                    self._data_outdir + outputvcf_filename,
                ],
                shell=False,
            )
        except Exception as e:
            print(e)
