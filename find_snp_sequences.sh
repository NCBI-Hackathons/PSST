#!/bin/bash

##################
#
# How to run this
#
#  for i in `cat snp_list.foo`; do echo $i; bash robot2.sh "$i"; done
#
##################

## it is worth noting that when selecting SRA accessions, SNPs that do not hit RNA are not useful for RNAseq or WXS.  Those SNPs can be filtered by this method
#  for i in `cat SNP_finder/snp_list.foo`; do echo $i; esearch -query "$i" -db snp | esummary | xtract -pattern DocumentSummary -element DOCSUM | egrep -o 'NM[^,]+' | sort -u; done 

## how to get accessions for RNA and protein from SNPs
#  esearch -query rs6003 -db snp | esummary | xtract -pattern DocumentSummary -element ACC

## how to get refseq amino acid changes from SNPs
#  esearch -query rs6003 -db snp | esummary | xtract -pattern DocumentSummary -element DOCSUM | egrep -o 'NP[^,]+' | sort -u

## Extract sequences from SNP

QUERY=$1
OUTPUT=$2

echo ">"$QUERY"_SEQ" >> $OUTPUT

esearch -query ${QUERY} -db snp | esummary | xtract -pattern DocumentSummary -element DOCSUM | egrep -o 'SEQ[^|]+' | sort -u | awk -F"SEQ=" '{print $2}' | sed 's/\[//' | sed 's/\/\w\]//' | sed 's/\///g' | sed 's/\]//g' | sed 's/\-//g'

#$OUTPUT | findFlankLengths.py >> flanking_lengths.txt

## For loop for extracting sequences from the SNP IDs, and compiling in fasta files

##for DISEASE in $(cat interesting_phenotypes);do SNP_IDs=$(echo ${DISEASE}_snps);FASTA_FILES=$(echo ${DISEASE}.fasta); for i in $(cat $SNP_IDs); do echo ">${i}_SEQ" >> $FASTA_FILES            
##> esearch -query $i -db snp | esummary | xtract -pattern DocumentSummary -element DOCSUM | egrep -o 'SEQ[^|]+' | sort -u | awk -F"SEQ=" '{print $2}' | sed 's/\[//' | sed 's/\/\w\]//' | sed 's/\///g' >> $FASTA_FILES


###

##makeblastdb -dbtype nucl -in SNP_SEQS.foo 


