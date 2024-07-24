#!/bin/bash
# this conversion only takes care of potential snp's (but not indels) reflected in the resulting fasta file.

refgenomefile=$1
inputvcffile=$2
resultfilepath=$3
base_dir=$4

#bcftools mpileup -f $refgenomefile $inputbamfile | bcftools call -c | bcftools filter -s LowQual -e '%QUAL<20' - | vcfutils.pl vcf2fq | seqtk seq -A -l 70 - | sed -r "s/(>)(.*)/\1$bam.consensus/g" > $resultfilepath
cat $refgenomefile | bcftools consensus $inputvcffile > $resultfilepath
