__author__ = 'spark'

import sys
import subprocess
import os
import gzip
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

import vcf
import re

__all__ = [
    "SeqAnnotator",
]

#IN_FILE_FORMAT='gb' # 'embl'
IN_FILE_FORMAT='embl'
OUT_FILE_FORMAT='embl'
EMBL_FORMAT='embl'
GB_FORMAT='gb'

class SeqAnnotator():
    _in_file_format = EMBL_FORMAT
    _out_file_format = EMBL_FORMAT
    # expects refseq_file in embl format and conseq_file in fasta format
    def __init__(self, refseq_file, conseq_file, data_outdir='data/out', ref_format=EMBL_FORMAT, out_format=EMBL_FORMAT):
        self._in_file_format = ref_format
        self._out_file_format = out_format
        # read the ref sequence file
        print(f'SeqAnnotator init: parsing the ref and consensus sequences: {refseq_file}, {conseq_file}')
        with open(refseq_file, 'r') as f:
            refseq_record = list(SeqIO.parse(f, self._in_file_format))[0]
        with open(conseq_file, 'r') as f:
            conseq_record = list(SeqIO.parse(f, "fasta"))[0]
        print('SeqAnnotator init: parsing done')

        self._refseq_file = refseq_file
        self._conseq_file = conseq_file
        self._conseq_file_nameonly = re.search('(.*)\\.(gb|embl|fasta)$', conseq_file).groups()[0]
        self._data_outdir = data_outdir

    #def executeRATT(self, prefix_str='consensus_ratt', transfer_type='Strain', ratt_home_str='/vagrant/Documents/workspace/GenomeAssemblyPipeline/lib/ratt'):
    def executeRATT(self, refseq_path=None, prefix_str='consensus_ratt', transfer_type='Strain', ratt_home_str='toolset/ratt', final_copyto_path=None):
        if (refseq_path == None):
            refseq_path = self._refseq_file
        pwd = os.path.realpath('.')
        os.environ['RATT_HOME'] = pwd + '/' + ratt_home_str
        os.environ['RATT_CONFIG'] = pwd + '/' + ratt_home_str + '/' + 'RATT.config'
        os.environ['RATT_VERBOSE'] = 'verbose'
        print(pwd)
        if (not os.path.exists(self._data_outdir)):
            os.makedirs(self._data_outdir)
        os.chdir(self._data_outdir)
        refseq_dir = os.path.dirname(refseq_path)
        print(pwd + '/' + refseq_dir)
        print(pwd + '/' + self._conseq_file)
        print(pwd + '/' + self._refseq_file)
        command = os.environ['RATT_HOME'] + '/' + 'start.ratt.sh'

        #import ipdb; ipdb.set_trace()
        subprocess.call([command, pwd + '/' + refseq_dir, pwd + '/' + self._conseq_file, str(prefix_str), transfer_type])
        #subprocess.call([command, pwd + '/' + self._refseq_file, pwd + '/' + self._conseq_file, prefix_str, transfer_type])

        if (final_copyto_path):
            subprocess.call('cp '+ pwd+'/'+self._data_outdir+'/consensus_ratt*.final.'+self._out_file_format + ' ' +  pwd + '/' + final_copyto_path +'/' + self._conseq_file_nameonly + '.' +self._out_file_format, shell=True)
        else:
            subprocess.call('cp '+ pwd+'/'+self._data_outdir+'/consensus_ratt*.final.'+self._out_file_format + ' ' +  pwd + '/' + self._data_outdir+'/' + self._conseq_file_nameonly + '.' +self._out_file_format, shell=True)
        os.chdir(pwd)

    def _isTrueVariantFeature(self, in_embl_seq, ref_embl_seq, var_record, depth_cutoff=20):
        start_seq_idx = var_record.POS - 1
        end_seq_idx = start_seq_idx + max(len(var_record.REF[0]),len(var_record.ALT[0]))

        ref_str = str(ref_embl_seq.seq[start_seq_idx:end_seq_idx])
        var_str = str(in_embl_seq.seq[start_seq_idx:end_seq_idx])

        if (var_str.lower() != ref_str.lower() and var_record.INFO['DP'] > depth_cutoff):
            return True
        else:
            return False

    def fillSeqGaps(self, in_embl_seq, ref_embl_seq):
        in_embl_str = str(in_embl_seq.seq) #.tostring()

        cur_idx = 0

        newfeatures = list(in_embl_seq.features) # copy the current feature list in the embl file
        print(f"fillSeqGaps: type: {str(type(in_embl_seq.seq))}")
        in_mutable_seq = MutableSeq(in_embl_seq.seq)

        rematch = re.search('[^ACGTacgt]+', in_embl_str)
        while (rematch != None):
            startidx = rematch.start()
            endidx = rematch.end()

            matchstr = str(rematch.group())

            if (cur_idx + startidx > 0):
                stridx1 = cur_idx + startidx -1
            else:
                stridx1 = cur_idx + startidx

            if (cur_idx + endidx < len(in_embl_str) - 1):
                stridx2 = cur_idx + endidx + 1
            else:
                stridx2 = cur_idx + endidx

            print(in_embl_str[stridx1:stridx2])
            seqidx1 = cur_idx+startidx
            seqidx2 = cur_idx+endidx
            in_mutable_seq[seqidx1:seqidx2] = ref_embl_seq.seq[seqidx1:seqidx2]

            # search for the coincident gene or CDS feature
            gene_feature = [f for f in in_embl_seq.features if ((seqidx1+1 in f or seqidx2+1 in f) and (f.type == 'gene' or f.type == 'CDS'))] 
            if (len(gene_feature) > 0):
                strandval = gene_feature[0].strand
            else:
                strandval = +1

            # make a new gap feature annotation
            gapfeature = SeqFeature(FeatureLocation(seqidx1, seqidx2, strand=strandval), type='gap', qualifiers={'original_gap_seq':matchstr, 'ref_seq':str(in_mutable_seq[seqidx1:seqidx2])})

            newfeatures.append(gapfeature)

            cur_idx = cur_idx + endidx
            rematch = re.search('[^ACGTacgt]+', in_embl_str[cur_idx:])

        #in_embl_seq.seq = in_mutable_seq.toseq() 
        # now populate the feature with the new feature list
        new_embl_seq = SeqRecord(in_mutable_seq, id=in_embl_seq.id, name=in_embl_seq.name, description=in_embl_seq.description)
        new_embl_seq.features = newfeatures

        print(f"fillSeqGaps: returning type: {str(type(new_embl_seq))}")

        return new_embl_seq

    def transferAnnotations(self, in_seq, ref_seq):
        ref_features = ref_seq.features
        newfeatures = list(in_seq.features)
        for f in ref_features:
            f_copy = f._shift(0)
            newfeatures.append(f_copy)

        # now populate the feature with the new feature list
        in_seq.features = newfeatures

        return in_seq

    def parseSeqRecordFile(self, seq_file):
        seq_format = re.search('(gb|embl|fasta)$', seq_file).group()
        with open(seq_file, 'r') as f:
            seq = list(SeqIO.parse(f, seq_format))[0]

        return seq


    def transferVCFAnnotations(self, in_seq, ref_seq, in_vcf_file):
        print(f"transferVCFAnnotations: reading a vcf file: {in_vcf_file}")
        with open(in_vcf_file, 'r') as vcff:
            vcf_reader = vcf.Reader(vcff)
            newfeatures = list(in_seq.features) # copy the current feature list in the embl file
            for var_record in vcf_reader:
                var_of_gene = [f for f in in_seq.features if (var_record.POS in f and (f.type == 'gene' or f.type == 'CDS'))] 
                var_of_rRNA = None
                if (len(var_of_gene) > 0):
                    var_of_rRNA = [f for f in in_seq.features if var_of_gene[0].location.start == f.location.start and var_of_gene[0].location.end == f.location.end and f.type == 'rRNA']

                isValidVariant = self._isTrueVariantFeature(in_seq, ref_seq, var_record)
                if (isValidVariant):
                    info_str = 'REF:' + str(var_record.REF) + ', ALT:' + str(var_record.ALT) + ', QAUL:' + str(var_record.QUAL) + ', Depth:' + str(var_record.INFO['DP'])
                else: # the variant is likely to be a false alarm
                    info_str = 'Invalidated variant call: REF:' + str(var_record.REF) + ', ALT:' + str(var_record.ALT) + ', QAUL:' + str(var_record.QUAL) + ', Depth:' + str(var_record.INFO['DP'])

                if (isValidVariant):
                    start_location = var_record.POS - 1
                    end_location = var_record.POS - 1 + len(var_record.ALT[0])

                    if (len(var_of_gene) > 0):
                        #gene_str = 'gene: ' + str(var_of_gene[0].qualifiers['gene'])
                        gene_str = 'label: ' + str(var_of_gene[0].qualifiers['gene'])
                        #print(f"transferVCFAnnotations: gene {gene_str} position: {var_of_gene[0].location.start}")
                        gene_si = var_of_gene[0].location.start
                        gene_ei = var_of_gene[0].location.end
                        if (var_of_gene[0].strand == 1):
                            rel_idx = start_location - gene_si
                            codon_idx = rel_idx%3
                            if (codon_idx == 0):
                                codon_si = start_location
                                codon_ei = start_location+3
                            elif (codon_idx == 1):
                                codon_si = start_location-1
                                codon_ei = start_location+2
                            elif (codon_idx == 2):
                                codon_si = start_location-2
                                codon_ei = start_location+1
                        else:
                            rel_idx = gene_ei -1 - start_location
                            codon_idx = rel_idx%3
                            if (codon_idx == 0):
                                codon_si = start_location-2
                                codon_ei = start_location+1
                            elif (codon_idx == 1):
                                codon_si = start_location-1
                                codon_ei = start_location+2
                            elif (codon_idx == 2):
                                codon_si = start_location
                                codon_ei = start_location+3

                        if (var_of_gene[0].strand == 1):
                            ref_codon_str = ref_seq.seq[codon_si:codon_ei]
                            alt_codon_str = in_seq.seq[codon_si:codon_ei]
                        else:
                            ref_codon_str = ref_seq.seq[codon_si:codon_ei].reverse_complement()
                            alt_codon_str = in_seq.seq[codon_si:codon_ei].reverse_complement()

                        ref_aa_str = ref_codon_str.translate()
                        alt_aa_str = alt_codon_str.translate()
                        mutation_note = ''
                        if (len(var_of_rRNA) == 0):
                            if (str(ref_aa_str) == str(alt_aa_str)):
                                mutation_note = 'silent mutation: ' + str(ref_aa_str) + ' -> ' + str(alt_aa_str)
                            else:
                                mutation_note = 'non-silent mutation: ' + str(ref_aa_str) + ' -> ' + str(alt_aa_str)

                        feature = SeqFeature(FeatureLocation(start_location, end_location, strand=+1), type=var_record.var_type, qualifiers={'gene':gene_str, 'info':info_str, 'note':mutation_note})
                        #feature.qualifiers['gene'] = var_of_gene
                    else:
                        feature = SeqFeature(FeatureLocation(start_location, end_location, strand=+1), type=var_record.var_type, qualifiers={'info':info_str})

                    newfeatures.append(feature)
                    # end of if
                # end of if

        # now populate the feature with the new feature list
        in_seq.features = newfeatures

        return in_seq

    def addHeaderAnnotations(self, embl_seq_records):
        newannotations = dict(embl_seq_records.annotations) # copy the current annotation list

        # add header annotations

        newannotations['comment']  = "CBCB's BSB1 strain of B. subtilis: sequenced using Illumina and assembled against the reference strain 168 (AL009126.3) dated 26-Feb-2014. The assembly was done using bowtie-2, and the varaint analysis was done using GATK (Genome Analysis ToolKit). The annotations (gene, CDS, tRNA and rRNA) from the reference strain AL009126 were transferred using RATT (Rapid Annotation Transfer Tool). The VCF file from the GATK pipeline was used to add annotations on variant calls."
        newannotations['organism'] = "Bacillus subtilis subsp. BSB1 from CBCB"
        newannotations['taxonomy'] = ['Bacteria', 'Firmicutes', 'Bacilli', 'Bacillales', 'Bacillaceae', 'Bacillus']
        newannotations['molecule_type'] = "DNA"


        source_feature = [f for f in embl_seq_records.features if (f.type == 'source')]
        if (len(source_feature) > 0):
            source_feature[0].qualifiers['db_xref']  = ['taxon:unknown']
            source_feature[0].qualifiers['organism'] = ['Bacillus subtilis subsp. BSB1']
            source_feature[0].qualifiers['strain'] = ['BSB1 (CBCB)']
            source_feature[0].qualifiers['sub_species'] = ['subtilis']


        embl_seq_records.annotations = newannotations
        embl_seq_records.description = 'Bacillus subtilis subsp. BSB1 (partially) complete genome.'

        return embl_seq_records

    def writeEMBLFile(self, embl_seq_records, out_embl_filename, out_format=None):
        # save the modified embl records
        SeqIO.write(embl_seq_records, self._data_outdir + '/' + out_embl_filename, self._out_file_format if out_format == None else out_format )



'''    def transferVCFAnnotation(self, in_embl_file, in_vcf_file):
        vcf_reader = vcf.Reader(open(in_vcf_file, 'r'))
        in_embl_seq = SeqIO.parse(open(in_embl_file, 'r'), self._in_file_format).next()

        newfeatures = list(in_embl_seq.features) # copy the current feature list in the embl file
        for var_record in vcf_reader:
            var_of_gene = [f.qualifiers['gene'] for f in in_embl_seq.features if var_record.POS in f and f.type == 'gene']
            var_of_cds = [f.qualifiers['gene'] for f in in_embl_seq.features if var_record.POS in f and f.type == 'CDS']
            info_str = 'REF:' + str(var_record.REF) + ', ALT:' + str(var_record.ALT) + ', QAUL:' + str(var_record.QUAL) + ', Depth:' + str(var_record.INFO['DP'])
            start_location = var_record.POS - 1
            end_location = var_record.POS - 1 + len(var_record.ALT[0])
            gene_str = ''
            if (len(var_of_gene) > 0 or len(var_of_cds) > 0):
                gene_str = 'gene_antt: ' + str(var_of_gene) + ', cds_antt: ' + str(var_of_cds)

            if len(gene_str) > 0:
                feature = SeqFeature(FeatureLocation(start_location, end_location), type=var_record.var_type, strand=+1, qualifiers={'gene':gene_str, 'info':info_str})
                #feature.qualifiers['gene'] = var_of_gene
            else:
                feature = SeqFeature(FeatureLocation(start_location, end_location), type=var_record.var_type, strand=+1, qualifiers={'info':info_str})

            newfeatures.append(feature)

        # now populate the feature with the new feature list
        in_embl_seq.features = newfeatures

        return in_embl_seq
'''

