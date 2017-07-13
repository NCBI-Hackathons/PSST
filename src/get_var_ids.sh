#!/bin/bash
# Copyright: NCBI 2017
# Author: Anmol Vohra

if [ "$#" -ne 1 ]; then
	echo "Description: given a disease phenotype, this script retrieves variants (currently only SNPs) related"
	echo "to the phenotype."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [phenotype]"
	exit 0
fi

PHENOTYPE=$1

# Find variant IDs

esearch -query ${PHENOTYPE} -db pubmed | elink -target snp | esummary | xtract -pattern DocumentSummary -element SNP_ID
