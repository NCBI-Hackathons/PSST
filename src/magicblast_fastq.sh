#!/bin/bash
# Copyright: NCBI 2017
# Author: Sean La

if [ "$#" -ne 4 ]; then
	echo "Description: Given a FASTQ file and a BLAST database, this script runs Magic-BLAST" 
	echo "             on each SRA dataset."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [FASTQ file] [BLAST DB name] [output dir] [threads]"
	exit 0
fi

# Retrieve the command line arguments
FASTQ=$1
DB_NAME=$2
OUTPUT_DIR=$3
THREADS=$4

# This prevents ambiguous splicing from occuring in Magic-BLAST
export MAPPER_NO_OVERLAPPED_HSP_MERGED=1
BASENAME=`basename "${FASTQ}"`
OUTPUT_FILE=${OUTPUT_DIR}/${BASENAME/.*/.mbo}
magicblast -query ${FASTQ} -infmt fastq -db ${DB_NAME} -out ${OUTPUT_FILE} -outfmt tabular -parse_deflines T -num_threads ${THREADS}
