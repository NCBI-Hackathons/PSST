#!/bin/sh
# Copyright: 2017
# Authors: Anmol Vohra, Ben Busby, and Sean La

set -e

if [ "$#" -ne 2 ]; then
	echo "Description: given a list of variant accessions, this script retrieves their respective flanking"
	echo "sequences and outputs into the given output file."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [variant accessions file] [output file]"
	exit 0
fi

ACCESSIONS=$1
OUTPUT=$2

# Make the output file blank
echo "" > ${OUTPUT}

# Find the flanking sequence for each variant
for ACC in $(cat ${ACCESSIONS}); do
	printf "${ACC}=" >> ${OUTPUT}
	esearch -query ${ACC} -db snp | esummary | xtract -pattern DocumentSummary -element DOCSUM | egrep -o 'SEQ[^|]+' | sort -u | awk -F"SEQ=" '{print $2}' >> ${OUTPUT} 
done
