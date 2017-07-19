#!/bin/bash
# Copyright: NCBI 2017
# Authors: Chipo Mashayomombe, Sean La

set -e

if [ "$#" -ne 2 ]; then
	echo "Description: Given a phenotype, outputs a file containing SRA accessions"
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [phenotype] [output file]"
	exit 0
fi

PHENOTYPE=$1
OUTPUT=$2

esearch -query ${PHENOTYPE} -db sra | esummary | xtract -pattern DocumentSummary -ACC @acc -block DocumentSummary -element "&ACC" | sed 's/.*\(SRR.\)/\1/p'|  sed 's/.*\(DRR.\)/\1/p'|  sed 's/.*\(ERR.\)/\1/p' > ${OUTPUT}
