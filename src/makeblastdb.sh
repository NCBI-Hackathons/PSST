#!/bin/bash
# Copyright: NCBI 2017
# Authors: Sean La & Anmol Vohra

set -e

if [ "$#" -ne 2 ]; then
	echo "Description: this script creates a BLAST database given a FASTA file and a directory."
	echo "The name of the BLAST database will be the basename of the input FASTA file sans extensions."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [FASTA file] [BLASTDB directory]"
	exit 0
fi

# Assign the command line arguments some readable names
FASTA=$1
DIR=$2

# Get the basename sans extensions of the FASTA file
FASTA_BASE=`basename "${FASTA}"`
FASTA_NAME=$( printf ${FASTA_BASE} | awk -F"." '{ print $1 }' )

# Make the BLAST database in the given directory
cd ${DIR}
makeblastdb -dbtype nucl -in ${FASTA} -out ${FASTA_NAME}
