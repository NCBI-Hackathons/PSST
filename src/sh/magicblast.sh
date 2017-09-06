#!/bin/bash

# Exit this script immediately when an error occurs
set -e

# The basename of this script
BASENAME=`basename "$0"`

description() {
    echo "Description: Given a list of NGS datasets and a BLAST database, this script runs Magic-BLAST" 
    echo "             on the NGS datasets and BLAST database."
}

usage() {
    echo "Usage: ${BASENAME} [-h help and usage] [-i input list file] [-d BLAST database name] [-o output directory]"
    echo "                     [-p MBO output paths list path] [-t threads for Magic-BLAST run]"
    echo "                     [-c max number of Magic-BLAST runs] [-f Input list contains only FASTQ paths (flag)]"
    echo "                     [-a Input list contains only FASTA paths (flag)]" 
    echo "Notes: if -f and -a are unset, assumes the input list file contains only SRA accessions."
}

# Initialize empty variables
IS_FASTQ=0
IS_FASTA=0

# Command line arguments
while getopts ":hi:d:o:p:t:c:fa" OPT; do
    case ${OPT} in
        h)  
            description
            usage
            exit 0
            ;;
        i)
            INPUT=${OPTARG}
            ;;
        d)
            DB_NAME=${OPTARG}
            ;;
        o)
            OUTPUT_DIR=${OPTARG}
            ;;
        p)
            PATHS_FILE=${OPTARG}
            ;;
        t)
            THREADS=${OPTARG}
            ;;
        c)
            MAX_PROCS=${OPTARG}
            ;;
        f)
            IS_FASTQ=1 
            ;;
        a)
            IS_FASTA=1
            ;;
        \?)
            echo "Invalid option: -${OPTARG}" >&2
            echo "${BASENAME} error; exiting."
            exit 1
    esac
done

OPTS_INCOMPLETE=0 

if [ -z "${INPUT}" ]; then
    echo "Error: please provide an input list file."
    OPTS_INCOMPLETE=1
fi
if [ -z "${DB_NAME}" ]; then
    echo "Error: please specify a BLAST database."
    OPTS_INCOMPLETE=1
fi
if [ -z "${OUTPUT_DIR}" ]; then
    echo "Error: please specify an output directory."
    OPTS_INCOMPLETE=1
fi
if [ -z "${PATHS_FILE}" ]; then
    echo "Error: please specify a path for the MBO output paths list file."
    OPTS_INCOMPLETE=1
fi
if [ -z "${THREADS}" ]; then
    echo "Error: please specify the number of threads to give to each Magic-BLAST run."
    OPTS_INCOMPLETE=1
fi
if [ -z "${MAX_PROCS}" ]; then
    echo "Error: please specify the maximum number of concurrent Magic-BLAST runs."
    OPTS_INCOMPLETE=1
fi
# Only one input file type can be specified.
if [ ${IS_FASTQ} -eq 1 ] && [ ${IS_FASTA} -eq 1 ]; then
    echo "Error: please specify one input file type."
    OPTS_INCOMPLETE=1
fi
if [ ${OPTS_INCOMPLETE} -eq 1 ]; then
    echo "${BASENAME} error; exiting."
    exit 1
fi

# This prevents ambiguous splicing from occuring in Magic-BLAST
export MAPPER_NO_OVERLAPPED_HSP_MERGED=0

# Empty the paths file
echo "" > ${PATHS_FILE}

for LINE in $(cat ${INPUT}); do
    # Input is a list of FASTQ files
    if [ ${IS_FASTQ} -eq 1 ]
    then
        FASTQ_BASENAME=`basename "$LINE"`
        OUTPUT_FILE=${OUTPUT_DIR}/${FASTQ_BASENAME}.mbo
        magicblast -query ${LINE} -infmt fastq -db ${DB_NAME} -outfmt tabular -parse_deflines T -num_threads ${THREADS} | awk -F'\t' 'FNR > 3 { if ($2 != "-") {print} }' > ${OUTPUT_FILE} & 
    # Input is a list of FASTA files
    elif [ ${IS_FASTA} -eq 1 ]
    then
        FASTA_BASENAME=`basename "$LINE"`
        OUTPUT_FILE=${OUTPUT_DIR}/${FASTA_BASENAME}.mbo
        magicblast -query ${LINE} -infmt fasta -db ${DB_NAME} -outfmt tabular -parse_deflines T -num_threads ${THREADS} | awk -F'\t' 'FNR > 3 { if ($2 != "-") {print} }' > ${OUTPUT_FILE} &
    # Input is a list of SRA accessions
    else
        OUTPUT_FILE=${OUTPUT_DIR}/${LINE}.mbo
        magicblast -sra ${LINE} -db ${DB_NAME} -outfmt tabular -parse_deflines T -num_threads ${THREADS} | awk -F'\t' 'FNR > 3 { if ($2 != "-") {print} }' > ${OUTPUT_FILE} &
    fi
    # Add the path to the current MBO output file to the list
    echo ${OUTPUT_FILE} >> ${PATHS_FILE}
    while [ $(jobs -r | wc -l) -ge ${MAX_PROCS} ]; do sleep 1; done
done
# Wait for all processes to finish before exiting
wait
