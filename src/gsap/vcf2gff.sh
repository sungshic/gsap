#!/bin/bash
# script to convert .vcf file to .gff3 format

#set variables that do not vary across detected polymorphisms: name of the program used, and name of the reference sequence

software_type=$(cat $1| grep "##source="|cut -b 10-)
reference_sequence=$(cat $1| grep "##reference"|cut -b 13-)

#remove comment lines in .vcf file and pipe to a while loop

sed /^#/d $1| while read line;

#iterate through lines that describe polymorphisms

do
#set variables that vary between detected polymorphisms

length_of_feature=$(echo $line|sed 's/.*LEN=\([^;]*\);.*/\1/')	
polymorphism_type=$(echo $line|sed 's/.*TYPE=\([^;]*\);.*/\1/')
start_pos=$(echo $line|awk '{FS="\t"; OFS="\t";print $2}')
end_pos=$(($start_pos+$length_of_feature-1))
ref_alt=$(echo $line|awk '{FS="\t";print $4,"mutates to",$5}')

#print each polymorphism out in .gff3 format
#the type (polymorphism_type) is the same as the text listed in the "TYPE" field in the .vcf file
#the start position is the same as in the .vcf file; the end position is the (start position + (length of feature - 1))
#the description field contains a text description like so: "(reference sequence) mutates to (mutation)

echo $line|awk '{FS="\t"; OFS="\t";print "'"$reference_sequence"'","'"$software_type"'","'"$polymorphism_type"'","'"$start_pos"'","'"$end_pos"'",$6,".",$3,"'"$ref_alt"'"}'


done
