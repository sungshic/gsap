"""
Implements glue functions for sequence annotation in genome DNA data.

Functions perform sequence annotation by combining annotations from
a reference genome and custom annotations from variant calls.

Author: Sunny Park

"""

import os
import re
import subprocess
from typing import TypedDict

import vcf
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

__all__ = [
    "SeqAnnotator",
]

# IN_FILE_FORMAT='gb' # 'embl'
IN_FILE_FORMAT = "embl"
OUT_FILE_FORMAT = "embl"
EMBL_FORMAT = "embl"
GB_FORMAT = "gb"


class HeaderAnnotationType(TypedDict, total=False):
    """defines the header_annotations parameter type of addHeaderAnnotations."""

    comment: str
    organism: str
    taxonomy: str
    molecule_type: str


class SeqAnnotator:
    """
    Offers methods for annotating a genome sequence.

    Following are the details of constructor parameters to initialise
    key input data for GSAP's sequence annotation operations

    :param refseq_file: the file path of a reference genome in FASTA format.
    :param conseq_file: the file path of a consensus sequence in FASTA format.
    :param data_outdir: the output directory where resulting files stored.
    :param ref_format: the reference sequence data format 'embl' or 'gb'
    :param out_format: the output sequence data format 'embl' or 'gb'
    """

    # expects refseq_file in embl format and conseq_file in fasta format
    def __init__(
        self,
        refseq_file: str,
        conseq_file: str,
        data_outdir: str = "data/out",
        ref_format: str = EMBL_FORMAT,
        out_format: str = EMBL_FORMAT,
    ):
        """Constructor method."""
        self._in_file_format = ref_format
        self._out_file_format = out_format
        # read the ref sequence file
        # with open(refseq_file) as f:
        #    refseq_record = next(SeqIO.parse(f, self._in_file_format))
        # print("SeqAnnotator init: parsing done")

        self._refseq_file = refseq_file
        self._conseq_file = conseq_file
        re_match = re.search("(.*)\\.(gb|embl|fasta)$", conseq_file)
        if type(re_match) is not re.Match:
            raise Exception("conseq_file has unsupported extension")
        self._conseq_file_nameonly = re_match.groups()[0]
        self._data_outdir = data_outdir

    def executeRATT(
        self,
        refseq_path: str,
        prefix_str: str = "consensus_ratt",
        transfer_type: str = "Strain",
        ratt_home_str: str = "toolset/ratt",
        final_copyto_path: str | None = None,
    ) -> None:
        """Transfers annotations from a reference genome to the assembled genome."""
        refseq_path = self._refseq_file
        pwd = os.path.realpath(".")
        os.environ["RATT_HOME"] = pwd + "/" + ratt_home_str
        os.environ["RATT_CONFIG"] = pwd + "/" + ratt_home_str + "/" + "RATT.config"
        os.environ["RATT_VERBOSE"] = "verbose"
        print(pwd)
        if not os.path.exists(self._data_outdir):
            os.makedirs(self._data_outdir)
        os.chdir(self._data_outdir)
        refseq_dir = os.path.dirname(refseq_path)
        print(pwd + "/" + refseq_dir)
        print(pwd + "/" + self._conseq_file)
        print(pwd + "/" + self._refseq_file)
        command = os.environ["RATT_HOME"] + "/" + "start.ratt.sh"

        subprocess.call(
            [
                command,
                pwd + "/" + refseq_dir,
                pwd + "/" + self._conseq_file,
                str(prefix_str),
                transfer_type,
            ]
        )

        if final_copyto_path:
            subprocess.call(
                [
                    "/usr/bin/cp",
                    self._data_outdir
                    + "/consensus_ratt*.final."
                    + self._out_file_format,
                    final_copyto_path
                    + "/"
                    + self._conseq_file_nameonly
                    + "."
                    + self._out_file_format,
                ],
                shell=False,
            )
        else:
            subprocess.call(
                [
                    "/usr/bin/cp",
                    self._data_outdir
                    + "/consensus_ratt*.final."
                    + self._out_file_format,
                    self._data_outdir
                    + "/"
                    + self._conseq_file_nameonly
                    + "."
                    + self._out_file_format,
                ],
                shell=False,
            )
        os.chdir(pwd)

    def _isTrueVariantFeature(
        self,
        in_embl_seq: SeqRecord,
        ref_embl_seq: SeqRecord,
        var_record: vcf.model._Record,
        depth_cutoff: int = 20,
    ) -> bool:
        """Method to filter variant calls based on sequence depth."""
        start_seq_idx = var_record.POS - 1
        end_seq_idx = start_seq_idx + max(
            len(var_record.REF[0]), len(var_record.ALT[0])
        )

        ref_str = str(ref_embl_seq.seq[start_seq_idx:end_seq_idx])
        var_str = str(in_embl_seq.seq[start_seq_idx:end_seq_idx])

        if var_str.lower() != ref_str.lower() and var_record.INFO["DP"] > depth_cutoff:
            return True
        else:
            return False

    def fillSeqGaps(self, in_embl_seq: SeqRecord, ref_embl_seq: SeqRecord) -> SeqRecord:
        """Fills unknown nucleotide characters with features annotating such gaps."""
        in_embl_str = str(in_embl_seq.seq)  # .tostring()

        cur_idx = 0

        newfeatures = list(
            in_embl_seq.features
        )  # copy the current feature list in the embl file
        print(f"fillSeqGaps: type: {type(in_embl_seq.seq)!s}")
        in_mutable_seq = MutableSeq(in_embl_seq.seq)

        rematch = re.search("[^ACGTacgt]+", in_embl_str)
        while type(rematch) is re.Match:
            startidx = rematch.start()
            endidx = rematch.end()

            matchstr = rematch.group()

            if cur_idx + startidx > 0:
                stridx1 = cur_idx + startidx - 1
            else:
                stridx1 = cur_idx + startidx

            if cur_idx + endidx < len(in_embl_str) - 1:
                stridx2 = cur_idx + endidx + 1
            else:
                stridx2 = cur_idx + endidx

            print(in_embl_str[stridx1:stridx2])
            seqidx1 = cur_idx + startidx
            seqidx2 = cur_idx + endidx
            in_mutable_seq[seqidx1:seqidx2] = ref_embl_seq.seq[seqidx1:seqidx2]

            # search for the coincident gene or CDS feature
            gene_feature = [
                f
                for f in in_embl_seq.features
                if (
                    (seqidx1 + 1 in f or seqidx2 + 1 in f)
                    and (f.type == "gene" or f.type == "CDS")
                )
            ]
            if len(gene_feature) > 0:
                strandval = gene_feature[0].location.strand
            else:
                strandval = +1

            # make a new gap feature annotation
            gapfeature = SeqFeature(
                FeatureLocation(seqidx1, seqidx2, strand=strandval),
                type="gap",
                qualifiers={
                    "original_gap_seq": matchstr,
                    "ref_seq": str(in_mutable_seq[seqidx1:seqidx2]),
                },
            )

            newfeatures.append(gapfeature)

            cur_idx = cur_idx + endidx
            rematch = re.search("[^ACGTacgt]+", in_embl_str[cur_idx:])

        # in_embl_seq.seq = in_mutable_seq.toseq()
        # now populate the feature with the new feature list
        new_embl_seq = SeqRecord(
            in_mutable_seq,
            id=in_embl_seq.id,
            name=in_embl_seq.name,
            description=in_embl_seq.description,
        )
        new_embl_seq.features = newfeatures

        print(f"fillSeqGaps: returning type: {type(new_embl_seq)!s}")

        return new_embl_seq

    def transferAnnotations(self, in_seq: SeqRecord, ref_seq: SeqRecord) -> SeqRecord:
        """Transfers annotations from a reference genome to the assembled genome."""
        ref_features = ref_seq.features
        newfeatures = list(in_seq.features)
        for f in ref_features:
            f_copy = f._shift(0)
            newfeatures.append(f_copy)

        # now populate the feature with the new feature list
        in_seq.features = newfeatures

        return in_seq

    def parseSeqRecordFile(self, seq_file: str) -> SeqRecord | None:
        """Parses a sequence file into a biopython SeqRecord object."""
        re_match = re.search("(gb|embl|fasta)$", seq_file)
        if re_match:
            seq_format = re_match.group()
            with open(seq_file) as f:
                seq = next(SeqIO.parse(f, seq_format))

            return seq
        return None

    def transferVCFAnnotations(
        self, in_seq: SeqRecord, ref_seq: SeqRecord, in_vcf_file: str
    ) -> SeqRecord | None:
        """Transfers annotations from a vcf file to the assembled genome."""
        print(f"transferVCFAnnotations: reading a vcf file: {in_vcf_file}")
        with open(in_vcf_file) as vcff:
            vcf_reader = vcf.Reader(vcff)
            newfeatures = list(
                in_seq.features
            )  # copy the current feature list in the embl file
            if type(vcf_reader) is not vcf.parser.Reader:
                return None
            for var_record in vcf_reader:
                var_of_gene = [
                    f
                    for f in in_seq.features
                    if (var_record.POS in f and (f.type == "gene" or f.type == "CDS"))
                ]
                var_of_rRNA = []
                if len(var_of_gene) > 0:
                    var_of_rRNA = [
                        f
                        for f in in_seq.features
                        if var_of_gene[0].location.start == f.location.start
                        and var_of_gene[0].location.end == f.location.end
                        and f.type == "rRNA"
                    ]

                isValidVariant = self._isTrueVariantFeature(in_seq, ref_seq, var_record)
                if isValidVariant:
                    info_str = (
                        "REF:"
                        + str(var_record.REF)
                        + ", ALT:"
                        + str(var_record.ALT)
                        + ", QAUL:"
                        + str(var_record.QUAL)
                        + ", Depth:"
                        + str(var_record.INFO["DP"])
                    )
                else:  # the variant is likely to be a false alarm
                    info_str = (
                        "Invalidated variant call: REF:"
                        + str(var_record.REF)
                        + ", ALT:"
                        + str(var_record.ALT)
                        + ", QAUL:"
                        + str(var_record.QUAL)
                        + ", Depth:"
                        + str(var_record.INFO["DP"])
                    )

                if isValidVariant:
                    start_location = var_record.POS - 1
                    end_location = var_record.POS - 1 + len(var_record.ALT[0])

                    if len(var_of_gene) > 0:
                        gene_str = "label: " + str(var_of_gene[0].qualifiers["gene"])
                        gene_si = var_of_gene[0].location.start
                        gene_ei = var_of_gene[0].location.end
                        if var_of_gene[0].location.strand == 1:
                            rel_idx = start_location - gene_si
                            codon_idx = rel_idx % 3
                            if codon_idx == 0:
                                codon_si = start_location
                                codon_ei = start_location + 3
                            elif codon_idx == 1:
                                codon_si = start_location - 1
                                codon_ei = start_location + 2
                            elif codon_idx == 2:
                                codon_si = start_location - 2
                                codon_ei = start_location + 1
                        else:
                            rel_idx = gene_ei - 1 - start_location
                            codon_idx = rel_idx % 3
                            if codon_idx == 0:
                                codon_si = start_location - 2
                                codon_ei = start_location + 1
                            elif codon_idx == 1:
                                codon_si = start_location - 1
                                codon_ei = start_location + 2
                            elif codon_idx == 2:
                                codon_si = start_location
                                codon_ei = start_location + 3

                        if var_of_gene[0].location.strand == 1:
                            ref_codon_str = ref_seq.seq[codon_si:codon_ei]
                            alt_codon_str = in_seq.seq[codon_si:codon_ei]
                        else:
                            ref_codon_str = ref_seq.seq[
                                codon_si:codon_ei
                            ].reverse_complement()
                            alt_codon_str = in_seq.seq[
                                codon_si:codon_ei
                            ].reverse_complement()

                        ref_aa_str = ref_codon_str.translate()
                        alt_aa_str = alt_codon_str.translate()
                        mutation_note = ""
                        if len(var_of_rRNA) == 0:
                            if str(ref_aa_str) == str(alt_aa_str):
                                mutation_note = (
                                    "silent mutation: "
                                    + str(ref_aa_str)
                                    + " -> "
                                    + str(alt_aa_str)
                                )
                            else:
                                mutation_note = (
                                    "non-silent mutation: "
                                    + str(ref_aa_str)
                                    + " -> "
                                    + str(alt_aa_str)
                                )

                        feature = SeqFeature(
                            FeatureLocation(start_location, end_location, strand=+1),
                            type=var_record.var_type,
                            qualifiers={
                                "gene": gene_str,
                                "info": info_str,
                                "note": mutation_note,
                            },
                        )
                        # feature.qualifiers['gene'] = var_of_gene
                    else:
                        feature = SeqFeature(
                            FeatureLocation(start_location, end_location, strand=+1),
                            type=var_record.var_type,
                            qualifiers={"info": info_str},
                        )

                    newfeatures.append(feature)
                    # end of if
                # end of if

        # now populate the feature with the new feature list
        in_seq.features = newfeatures

        return in_seq

    def addHeaderAnnotations(
        self, embl_seq_records: SeqRecord, header_annotations: HeaderAnnotationType
    ) -> SeqRecord:
        """Adds header annotations to the final assembled, annotated sequence file."""
        newannotations = dict(
            embl_seq_records.annotations
        )  # copy the current annotation list

        # add header annotations

        newannotations["comment"] = header_annotations["comment"]
        newannotations["organism"] = header_annotations["organism"]
        newannotations["taxonomy"] = header_annotations["taxonomy"]
        newannotations["molecule_type"] = header_annotations["molecule_type"]

        embl_seq_records.annotations = newannotations
        embl_seq_records.description = (
            "Bacillus subtilis subsp. BSB1 (partially) complete genome."
        )

        return embl_seq_records

    def writeEMBLFile(
        self,
        embl_seq_records: SeqRecord,
        out_embl_filename: str,
        out_format: str | None = None,
    ) -> None:
        """Saves a biopython SeqRecord into an EMBL file."""
        SeqIO.write(
            embl_seq_records,
            self._data_outdir + "/" + out_embl_filename,
            self._out_file_format if out_format is None else out_format,
        )
