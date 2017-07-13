#!/bin/bash
# Copyright: NCBI 2017
# Author: Sean La

if [ "$#" -ne 4 ]; then
	echo "Description: given a file containing SRA accessions and a BLAST database, this script runs Magic-BLAST" 
	echo "on each SRA accessions."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [SRA accessions file] [BLASTDB dir path] [BLAST database name] [output directory]"
	exit 0
fi

# Maximum number of concurrent Magic-BLAST runs at a time
MAX_PROCS=10

# Retrieve the command line arguments
SRA=$1
BLASTDB=$2
DB_NAME=$3
OUTPUT_DIR=$4
MAPPER_NO_OVERLAPPED_HSP_MERGED=1

export BLASTDB
# This prevents ambiguous splicing from occuring in Magic-BLAST
export MAPPER_NO_OVERLAPPED_HSP_MERGED

for ACC in $(cat ${SRA}); do
	OUTPUT_FILE=${OUTPUT_DIR}/${ACC}.mbo
	magicblast -sra ${ACC} -db ${DB_NAME} -outfmt tabular -out ${OUTPUT_FILE} &
	# Limit the number of child processes running so we don't overload the local computer
	while [ $(jobs | wc -l) -gt "${MAX_PROCS}" ]; do sleep 1; done
done

# Wait for all processes to finish before exiting
wait
