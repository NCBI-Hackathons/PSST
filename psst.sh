#!/bin/bash
# Copyright: NCBI 2017
# Authors: Sean La

# Exit the script immediately if there is an error
set -e

if [ "$#" -ne 5 ]; then
	printf "Description: Given lists of SRA and SNP accessions, determines the set of SNPs that occur in each SRA"
	printf " dataset.\n"
	BASENAME=`basename "$0"`
	printf "Usage: ${BASENAME} [SRA accessions] [SNP accessions] [working directory] [threads] "
	printf "[max number of child processes]\n"
	exit 0
fi

## Retrieve the command line arguments and set up directories, paths
SRA_ACC=$1
SNP_ACC=$2
DIR=$3
THREADS=$4
PROCS=$5

mkdir -p ${DIR} # If the working directory does not exist, create it
SRC=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/src
export BLASTDB=${DIR}

## Find the variant flanking sequences
echo "Finding SNP flanking sequences..."
SNP_FLANKS=${DIR}/snp_flanks.txt
${SRC}/get_var_flanks.sh ${SNP_ACC} ${SNP_FLANKS}

## Get the variant information, i.e. the start and stop positions of the major allele variant in the flanking sequence
## and the length of the major allele
echo "Getting SNP flank information..."
SNP_INFO=${DIR}/snp_info.txt
${SRC}/find_var_info.py -i ${SNP_FLANKS} -o ${SNP_INFO} 

## Construct a FASTA file out of the variant flanking sequence file
echo "Constructing FASTA file out of the SNP flanking sequences..."
SNP_FASTA=${DIR}/snp_flanks.fasta
${SRC}/var_flanks_to_fasta.py -i ${SNP_FLANKS} -o ${SNP_FASTA}

## Create a BLAST database out of the variant FASTA file
echo "Creating a BLAST database out of the SNP flanks FASTA file..."
${SRC}/makeblastdb.sh ${SNP_FASTA} ${DIR}

## Align the SRA datasets onto the variants (a la the BLAST database) using Magic-BLAST
echo "Aligning SRA datasets onto the SNPs..."
MBO_DIR=${DIR}/mbo # We will store the .mbo files here
mkdir -p ${MBO_DIR} # Create the directory if it doesn't exist yet
${SRC}/magicblast.sh ${SRA_ACC} snps_flanks ${MBO_DIR} ${THREADS} ${PROCS}

## Call variants in the SRA datasets
echo "Calling SNPs..."
TSV=${DIR}/${PHENOTYPE}.tsv
${SRC}/call_variants.py -m ${MBO_DIR} -v ${SNP_INFO} -f ${SNP_FASTA} -o ${TSV}
echo "PSST complete. Result file can be found at:"
echo ${TSV}
