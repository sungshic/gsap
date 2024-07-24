__author__ = 'spark'

import os
import sh
import subprocess
from gsap import PreProcessor, VariantCaller, SeqAnnotator

# make a copy of the environment
env = dict(os.environ)
#env['CLASSPATH']
#subprocess.call(['java', '-jar', 'GATK/GenomeAnalysisTK.jar'])

base_dir = './'

print('starting the assembly pipeline')


pre_processor = PreProcessor(base_dir+'src/gsap/', 'AL009126_v11',
                             base_dir+'tests/data/sample1/SP2_S6_L001_R1_001.fastq.gz',
                             base_dir+'tests/data/sample1/SP2_S6_L001_R2_001.fastq.gz',
                             base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.fasta',
                             base_dir+'tests/data/test_ref_data/',
                             base_dir+'tests/data/tmp/')


variant_caller = VariantCaller(base_dir+'src/gsap/', base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.fasta',
                               base_dir+'tests/data/test_ref_data/addrg_reads.bam',
                               base_dir+'tests/data/test_ref_data/',
                               base_dir+'tests/data/tmp/')


#print('######## pre-processing raw reads... ########')
## run the pre-processing pipeline
#
#print('quality-based trimming')
#pre_processor.trimSeqReads()
#
#print('###############assembly de novo...')
#pre_processor.assembleDeNovo()
#
##print 'find a contig with EVO_insert'
##ins_blaster = InsertBlaster('_query.fasta', outdir_base+'denovo/spades/output/contigs.fasta')
##ins_blaster.executeBlast()
##ins_blaster.makeAnnotatedSeq('_psg1729_evot2_payload.embl', '_evo_insert_denovo.fasta', 'genomeassemblypipeline/data/gandalf/basespace/bsb1_evot2/assembly/')
#
#print('aligning PE reads to reference...')
#pre_processor.alignPairedEndReads('aligned_reads.sam')
#
#
#print('sorting aligned reads...')
#pre_processor.sortSAMIntoBAM(base_dir+'tests/data/test_ref_data/aligned_reads.sam', 'sorted_reads.bam') # this command requires lots of memory, 4096MB used
#
#print('build bam index for sorted_reads.bam')
#pre_processor.buildBamIndex(base_dir+'tests/data/test_ref_data/sorted_reads.bam', base_dir+'tests/data/test_ref_data/sorted_reads.bai')
#
#print('removing unmapped reads...')
#pre_processor.removeUnmappedReads(base_dir+'tests/data/test_ref_data/sorted_reads.bam', base_dir+'tests/data/test_ref_data/filtered_sample_aligned.bam', base_dir+'tests/data/test_ref_data/unmapped_reads.bam')
#
#print('build bam index for filtered_sample_aligned.bam')
#pre_processor.buildBamIndex(base_dir+'tests/data/test_ref_data/filtered_sample_aligned.bam', base_dir+'tests/data/test_ref_data/filtered_sample_aligned.bai')
#
#print('sorting the filtered bam...')
#pre_processor.sortBAMIntoBAM('tests/data/test_ref_data/filtered_sample_aligned.bam', 'filtered_sample_aligned_sorted.bam')
#print('marking duplicate reads...')
#pre_processor.markDuplicates(base_dir+'tests/data/test_ref_data/filtered_sample_aligned_sorted.bam', 'dedup_reads.bam', 'metrics.txt')
#print('adding read group...')
#pre_processor.addReadGroup(base_dir+'tests/data/test_ref_data/dedup_reads.bam', 'addrg_reads.bam')
#print('building bam index from read groups...')
#pre_processor.buildBamIndex(base_dir+'tests/data/test_ref_data/addrg_reads.bam', base_dir+'tests/data/test_ref_data/addrg_reads.bai')
#print('building reference fasta index...')
#pre_processor.buildRefFastaIndex(base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.fasta',
#                                 base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.dict')
#
#
#
#
#print('removing unmapped reads...')
##pre_processor.removeUnmappedReads('data/Illumina/Sample1BSB1/sample1_aligned.bam', 'data_output/clc/filtered_sample_aligned.bam')
#
#print('sorting the filtered bam...')
##pre_processor.sortBAMIntoBAM('data_output/filtered_sample_aligned.bam', 'sorted_reads.bam')
#print('marking duplicates from the filtered bam...')
##pre_processor.markDuplicates('data_output/sorted_reads.bam', 'dedup_reads.bam', 'metrics.txt')
#
#print('adding read group...')
##pre_processor.addReadGroup('data_output/dedup_reads.bam',
##                           'addrg_reads.bam')
#print('building bam index from read groups...')
##pre_processor.buildBamIndex('data_output/addrg_reads.bam')
#print('building reference fasta index...')
##pre_processor.buildRefFastaIndex('data/RefGenome/AL009126.fasta',
##                                 'data/RefGenome/AL009126.dict')
#print('calling variants...')
#variant_caller.callVariants()
#
#
## use glimmer to find putative ORFs
##subprocess.call(['tigr-glimmer', 'long-orfs', '-n', '-t' '1.15', 'data_output/clc/consensus.fasta', 'data_output/clc/consensus.longorfs'])
##subprocess.call(['./g3-from-scratch.sh', 'data_output/clc/consensus.fasta', 'data_output/clc/consensus_g3'])
##subprocess.call(['./g3-from-scratch.sh', 'data/RefGenome/AL009126.fasta', 'data_output/clc/refseq_g3'])


#print('finding a consensus fasta...')
#subprocess.call(['./src/gsap/getconsensusfasta3.sh', base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.fasta', base_dir+'tests/data/test_ref_data/raw_variants.vcf.gz', base_dir+'tests/data/test_ref_data/consensus.fasta', base_dir+'src/gsap'])
# fastq to fasta
# the fastq file need to be truncated using sed command 'sed -n '1,70262 p' addrg_reads_clc_cons.fq > consensus.fasta'
# the header needs to change from @... to >...

refseq_file = base_dir+'tests/data/RefGenome/AL009126_v11/AL009126.gb'
refseq_format = 'gb'
conseq_file = base_dir+'tests/data/test_ref_data/consensus.fasta'
outdir_base = base_dir+'tests/data/test_ref_data/'
seqannotator = SeqAnnotator(refseq_file, conseq_file, outdir_base, ref_format=refseq_format)

in_seq = seqannotator.parseSeqRecordFile(conseq_file)
ref_seq = seqannotator.parseSeqRecordFile(refseq_file)
print(f'transfering annotations... in_seq_type: {str(type(in_seq))} in_seq: {in_seq}')
in_seq_records = seqannotator.transferAnnotations(in_seq, ref_seq)
in_seq_records = seqannotator.fillSeqGaps(in_seq_records, ref_seq)
in_vcf_gz_file = outdir_base+'raw_variants.vcf.gz'
in_vcf_file = outdir_base+'raw_variants.vcf'
print(f'sh.gunzip: {in_vcf_gz_file}')
sh.gunzip('-kf', in_vcf_gz_file)
print(f'sh.gunzip: returned')
in_seq_records = seqannotator.transferVCFAnnotations(in_seq_records, ref_seq, in_vcf_file)
print(f'transferVCFAnnotations: returned')
in_seq_records = seqannotator.addHeaderAnnotations(in_seq_records)
print('writing genbank file...')

out_gb_file = 'BSB1_annotated.gb'
seqannotator.writeEMBLFile(in_seq_records, out_gb_file, 'gb')

'''
refseq_file = 'data/RefGenome/AL009126.embl'
conseq_file = 'data_output/168consensus.fasta'
data_outdir = 'data_output/dump'
seqannotator = SeqAnnotator(refseq_file, conseq_file, data_outdir)
seqannotator.executeRATT(ratt_home_str='/home/spark/workspace/genomeassemblypipeline/lib/ratt')

in_embl_file = 'data/RefGenome/AL009126.embl'
in_vcf_file = 'data_output/raw_variangts.vcf'
out_embl_file = '168_annotated.embl'

print 'transfering vcf annotations...'
embl_seq_records = seqannotator.transferVCFAnnotation(in_embl_file, in_vcf_file)
seqannotator.addHeaderAnnotations(embl_seq_records)
print 'writing embl file...'
seqannotator.writeEMBLFile(embl_seq_records, out_embl_file)
'''

print('finished assembling the genome sequence')
