#!/bin/bash

set -e

if [ $# -ne 3 ]; then
	echo "Desciption: Calculates expression of RNA-seq data"
	BASENAME=`basename $0`
	echo "Usage: ${BASENAME} [RNA-seq FASTQ file] [Exome FASTA file] [Working dir] [Output path]"
	exit 0
fi

RNA=$1
EXOME=$2
DIR=$3
OUTPUT=$4

export BLASTDB=${DIR}

## Remove Poly(A) and Poly(T) tails from RNA reads
TAILESS=${DIR}/rna-no_tails-good.fastq
BAD=${DIR}/rna-no_tails-bad.fastq
prinseq-lite.pl -fastq ${RNA} -out_format 1 -trim_tail_left 5 -trim_tail_right 5 -out_good ${TAILESS} -out_bad ${BAD} 

## Mask low complexity regions in the exome BLAST database 
MASKED_EXOME=${DIR}/masked_exome.fasta
dustmasker -in ${EXOME} -outfmt fasta -out ${MASKED_EXOME} 

## Create a blast database out of the masked exome
${PWD}/makeblastdb.sh ${MASKED_EXOME} ${DIR}

## Align the RNA-seq data onto the exome
ALIGNED_RNA=${DIR}/aligned_rna_onto_exome.mbo
magicblast -query ${EXOME} -db masked_exome -outfmt tabular -out ${ALIGNED_RNA} -lcase_masking -splice T

## Calculate expression
${PWD}/calculate_expression.py -i ${ALIGNED_RNA} -o ${OUTPUT}

echo "RNA expression calculation complete."
