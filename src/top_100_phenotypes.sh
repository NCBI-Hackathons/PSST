#!/bin/bash
# Copyright: NCBI 2017
# Authors: Chipo Mashayomombe, Sean La

if [ "$#" -ne 1 ]; then
	echo "Description: Finds the top 100 phenotypes described in the latest release of ClinVar for human"
	echo "             genome assembly GRCh38."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [working directory]"
	exit 0
fi

# Retrieve the command line arguments
DIR=$1

if [ ! -f ${DIR}/top_100_disease_phenotypes.txt ]; then
	# Get the latest GRCh38 release of ClinVar
	if [ ! -f ${DIR}/GRCh38_latest_clinvar.vcf ]; then
		wget -P ${DIR} ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
		gzip -d ${DIR}/clinvar.vcf.gz
		mv ${DIR}/clinvar.vcf ${DIR}/GRCh38_latest_clinvar.vcf
	fi
	# Get the top 100 disease phenotypes 
	awk '{if ($0 ~ /CLNSIG=5/) print $0}' ${DIR}/GRCh38_latest_clinvar.vcf| egrep -o 'CLNDBN[^:]+' | awk -F";" '{print $1}' | sort | uniq -c | sort -nr | grep -v "not_provided" | tr '|' '\n' | tr '\' '\n' | sed 's/CLNDBN[=]//g' | sort -nr  | head -100 > ${DIR}/top100_disease_phenotypes.txt
fi
