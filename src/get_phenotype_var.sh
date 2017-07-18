#!/bin/bash
# Copyright: NCBI 2017
# Author: Anmol Vohra

set -e

if [ "$#" -ne 2 ]; then
	echo "Description: given a disease phenotype, this script retrieves variants (currently only SNPs) related"
	echo "to the phenotype."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [phenotype] [output file]"
	exit 0
fi

PHENOTYPE=$1
OUTPUT=$2

# Find variant IDs

esearch -query ${PHENOTYPE} -db pubmed | elink -target snp | esummary | xtract -pattern DocumentSummary -element SNP_ID > ${OUTPUT}
