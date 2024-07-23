#!/bin/bash
# this conversion only takes care of potential snp's (but not indels) reflected in the resulting fasta file.

refgenomefile=$1
inputbamfile=$2
resultfilepath=$3
base_dir=$4

samtools mpileup -uf $refgenomefile $inputbamfile | bcftools call -c | bcftools filter -s LowQual -e '%QUAL<20' - | $4/toolset/bcftools/misc/vcfutils.pl vcf2fq | seqtk seq -A -l 70 - | sed -r "s/(>)(.*)/\1$bam.consensus/g" > $resultfilepath
