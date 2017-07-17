#!/bin/bash
# Copyright: NCBI 2017
# Author: Sean La

if [ "$#" -ne 3 ]; then
	echo "Description: Given a file containing SRA accessions and a BLAST database, this script runs Magic-BLAST" 
	echo "             on each SRA dataset."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [SRA accessions file] [BLAST database name] [output directory]"
	exit 0
fi

# Maximum number of concurrent Magic-BLAST runs at a time
MAX_PROCS=10

# Retrieve the command line arguments
SRA=$1
DB_NAME=$2
OUTPUT_DIR=$3

# This prevents ambiguous splicing from occuring in Magic-BLAST
export MAPPER_NO_OVERLAPPED_HSP_MERGED=1

for ACC in $(cat ${SRA}); do
	OUTPUT_FILE=${OUTPUT_DIR}/${ACC}.mbo
	magicblast -sra ${ACC} -db ${DB_NAME} -outfmt tabular -out ${OUTPUT_FILE} &
	# Limit the number of child processes running so we don't overload the local computer
	while [ $(jobs | wc -l) -gt "${MAX_PROCS}" ]; do sleep 1; done
done

# Wait for all processes to finish before exiting
wait
