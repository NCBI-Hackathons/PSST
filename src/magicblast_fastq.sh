#!/bin/bash
# Copyright: NCBI 2017
# Author: Sean La

if [ "$#" -ne 4 ]; then
	echo "Description: Given a FASTQ file and a BLAST database, this script runs Magic-BLAST" 
	echo "             on each SRA dataset. Assumes the BLASTDB exists in the output directory."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [FASTQ file] [BLAST DB name] [output dir] [threads]"
	exit 0
fi

# Retrieve the command line arguments
FASTQ=$1
DB_NAME=$2
OUTPUT_DIR=$3
THREADS=$4

# Set the BLASTDB path to the working directory
export BLASTDB=${OUTPUT_DIR}

# This prevents ambiguous splicing from occuring in Magic-BLAST
export MAPPER_NO_OVERLAPPED_HSP_MERGED=1

MBO_DIR=${OUTPUT_DIR}/mbo

BASENAME=`basename "${FASTQ}"`
OUTPUT_FILE=${MBO_DIR}/${BASENAME/.*/.mbo}
magicblast -query ${FASTQ} -infmt fastq -db ${DB_NAME} -out ${OUTPUT_FILE} -outfmt tabular -parse_deflines T -num_threads ${THREADS}
