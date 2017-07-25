#!/bin/bash
# Copyright: NCBI 2017
# Authors: Sean La

# Exit the script immediately if there is an error
set -e

if [ "$#" -ne 2 ]; then
	echo "Description: Given a disease/phenotype, this script retrieves all related variant and SRA accession"
	echo "             numbers and determines the set of variants that occur in each SRA dataset."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [disease] [working directory] [threads per Magic-BLAST run] [max num of child procs]"
	exit 0
fi

## Retrieve the command line arguments and set up directories, paths
PHENOTYPE=$1
DIR=$2
THREADS=$3
PROCS=$4

mkdir -p ${DIR} # If the working directory does not exist, create it
SRC=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/src
export BLASTDB=${DIR}

## Find the variant and SRA accessions
echo "Finding variant and SRA accessions..."
VAR_ACC=${DIR}/var_acc.txt
${SRC}/get_phenotype_var.sh ${PHENOTYPE} ${VAR_ACC} 
SRA_ACC=${DIR}/sra_acc.txt
${SRC}/get_phenotype_sra.sh ${PHENOTYPE} ${SRA_ACC}

## Find the variant flanking sequences
echo "Finding variant flanking sequences..."
VAR_FLANKS=${DIR}/var_flanks.txt
${SRC}/get_var_flanks.sh ${VAR_ACC} ${VAR_FLANKS}

## Get the variant information, i.e. the start and stop positions of the major allele variant in the flanking sequence
## and the length of the major allele
echo "Getting variant information..."
VAR_INFO=${DIR}/var_info.txt
${SRC}/find_var_info.py -i ${VAR_FLANKS} -o ${VAR_INFO} 

## Construct a FASTA file out of the variant flanking sequence file
echo "Constructing FASTA file out of the variant flanking sequences..."
VAR_FASTA=${DIR}/variants.fasta
${SRC}/var_flanks_to_fasta.py -i ${VAR_FLANKS} -o ${VAR_FASTA}

## Create a BLAST database out of the variant FASTA file
echo "Creating a BLAST database out of the variant FASTA file..."
${SRC}/makeblastdb.sh ${VAR_FASTA} ${DIR}

## Align the SRA datasets onto the variants (a la the BLAST database) using Magic-BLAST
echo "Aligning SRA datasets onto the variants..."
MBO_DIR=${DIR}/mbo # We will store the .mbo files here
mkdir -p ${MBO_DIR} # Create the directory if it doesn't exist yet
${SRC}/magicblast.sh ${SRA_ACC} variants ${MBO_DIR} ${THREADS} ${PROCS}

## Call variants in the SRA datasets
echo "Calling variants..."
TSV=${DIR}/${PHENOTYPE}.tsv
${SRC}/call_variants.py -m ${MBO_DIR} -v ${VAR_INFO} -f ${VAR_FASTA} -o ${TSV}
echo "PSST pipeline done. Result file can be found at ${TSV}."
