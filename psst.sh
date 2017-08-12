#!/bin/bash
# Copyright: NCBI 2017
# Authors: Sean La

# Exit the script immediately if there is an error
set -e

# Description and usage messages
description() {
    printf "PSST version 1.1\n"
	printf "Description: Given lists of SRA and SNP accessions, determines the set of SNPs that occur in each SRA"
	printf " dataset.\n"
}

usage() { 
	BASENAME=`basename "$0"`
	printf "Usage: ${BASENAME} [-h description and usage] [-s SRA accessions] [-n SNP accessions]\n"
    printf "               [-f FASTQ file] [-d working directory] [-e email for Entrez]\n"
    printf "               [-t threads] [-p max number of child processes]\n"
    echo "Notes:"
    echo "'-h' is the only non-mandatory parameter."
    echo "Exactly one of '-s' or '-f' must be provided as an argument."
    echo "All other arguments are mandatory."
}

# Command line arguments
while getopts ":hs:f:n:d:e:t:p:" opt; do
    case ${opt} in
        h)
            description 
            usage
            exit 0
            ;;
        s) # path to the SRA accessions file
            SRA_ACC=${OPTARG}
            ;;
        f) # path to FASTQ file
            FASTQ=${OPTARG}
            ;;
        n) # path to the SNP accessions file
            SNP_ACC=${OPTARG}
            ;;
        d) # path to the working directory
            DIR=${OPTARG}
            ;;
        e) # email address to give to Entrez
            EMAIL=${OPTARG}
            ;;
        t) # number of threads per Magic-BLAST run
            THREADS=${OPTARG}
            ;;
        p) # maximum number of child processes
            PROCS=${OPTARG}
            ;;
        \?)
            echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
    esac
done

# Check whether all the necessary command line options have all been inputted and are correct
if [ -z "${SRA_ACC}" ] && [ -z "${FASTQ}" ]; then
    echo "Error: please provide either an SRA accessions file or a FASTQ file."
    OPTS_INCOMPLETE=0
fi
if [ ! -z "${SRA_ACC}" ] && [ ! -z "${FASTQ}" ]; then
    echo "Error: please provide only one of either an SRA accessions file or a FASTQ file."
    OPTS_INCOMPLETE=0
fi
if [ -z "${SNP_ACC}" ]; then
    echo "Error: please provide a SNP accessions file."
    OPTS_INCOMPLETE=0
fi
if [ -z "${DIR}" ]; then
    echo "Error: please specify a working directory."
    OPTS_INCOMPLETE=0
fi
if [ -z "${EMAIL}" ]; then
    echo "Error: please provide an email address for Entrez."
    OPTS_INCOMPLETE=0
fi
if [ -z "${THREADS}" ]; then
    echo "Error: please specify the number of threads to give to each Magic-BLAST run."
    OPTS_INCOMPLETE=0
fi
if [ -z "${PROCS}" ]; then
    echo "Error: please specify the maximum number of child processes for this program."
    OPTS_INCOMPLETE=0
fi
# Exit the script if the command line options are incomplete or incorrect
if [ -n "${OPTS_INCOMPLETE}" ]; then
    exit 1
fi

## Retrieve the command line arguments and set up directories, paths
mkdir -p ${DIR} # If the working directory does not exist, create it
SRC=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/src
export BLASTDB=${DIR}

## Find the variant flanking sequences
echo "Finding SNP flanking sequences..."
SNP_FLANKS=${DIR}/snp_flanks.txt
${SRC}/get_var_flanks.py -i ${SNP_ACC} -e ${EMAIL} -o ${SNP_FLANKS}

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

# Either run Magic-BLAST on list of SRA accessions or on the single FASTQ file
if [ -n "${SRA_ACC}" ]; then
    ${SRC}/magicblast_sra.sh ${SRA_ACC} snp_flanks ${MBO_DIR} ${THREADS} ${PROCS}
else
    ${SRC}/magicblast_fastq.sh ${FASTQ} snp_flanks ${MBO_DIR} ${THREADS}
fi

## Call variants in the SRA datasets
echo "Calling SNPs..."
TSV=${DIR}/results.tsv
declare -i COMBINED_PROCS
COMBINED_PROCS=${THREADS}*${PROCS}
${SRC}/call_variants.py -m ${MBO_DIR} -v ${SNP_INFO} -f ${SNP_FASTA} -p ${COMBINED_PROCS} -o ${TSV}
echo "PSST run complete. Result file can be found at:"
echo ${TSV}
