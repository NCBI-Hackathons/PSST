#!/bin/bash
# Copyright: NCBI 2017
# Author: Sean La

if [ "$#" -ne 5 ]; then
	echo "Description: Given a file containing SRA accessions and a BLAST database, this script runs Magic-BLAST" 
	echo "             on each SRA dataset. Assumes the BLASTDB is in the output directory."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [SRA accessions file] [BLAST DB name] [output dir] [threads] [max child procs]"
	exit 0
fi

# Retrieve the command line arguments
SRA=$1
DB_NAME=$2
OUTPUT_DIR=$3
THREADS=$4
MAX_PROCS=$5

# Sets the BLASTDB path to the output directory
export BLASTDB=${OUTPUT_DIR}

# This prevents ambiguous splicing from occuring in Magic-BLAST
export MAPPER_NO_OVERLAPPED_HSP_MERGED=1

# The directory to store the mbo files
MBO_DIR=${OUTPUT_DIR}/mbo

for ACC in $(cat ${SRA}); do
	OUTPUT_FILE=${MBO_DIR}/${ACC}.mbo
	magicblast -sra ${ACC} -db ${DB_NAME} -outfmt tabular -parse_deflines T -num_threads ${THREADS} | awk -F'\t' 'FNR > 3 { if ($2 != "-") { print $2,$9,$10,$17 } }' > ${OUTPUT_FILE} &
	# Limit the number of child processes running so we don't overload the local computer
	while [ $(jobs | wc -l) -ge "${MAX_PROCS}" ]; do sleep 1; done
done

# Wait for all processes to finish before exiting
wait
