#!/bin/bash
# Copyright: NCBI 2017
# Author: Chipo Mashayomombe, Sean La

set -e

if [ "$#" -ne 1 ]; then
	echo "Description retrieves the latest ClinVar GRCh38 VCF file release, takes the top 100 phenotypes"
	echo "depicted in the VCF file and outputs them into another file."
	BASENAME=`basename "$0"`
	echo "${BASENAME} [working directory]"
	exit 0
fi

DIR=$1
VCF=${DIR}/GRCh38_latest_clinvar.vcf
TOP_100=${DIR}/top100_disease_phenotypes.txt
TEMP_FILE=${DIR}/temp.txt

# Retrieve the latest ClinVar release for human genome assembly GRCh38
if [ ! -f ${VCF} ]; then 
	wget -P ${DIR} ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
	gzip -d ${DIR}/clinvar.vcf.gz
	mv ${DIR}/clinvar.vcf ${VCF}
fi

awk '{if ($0 ~ /CLNSIG=5/) print $0}' ${VCF} | egrep -o 'CLNDBN[^:]+' | awk -F";" '{print $1}' | sort | uniq -c | sort -nr | grep -v "not_provided" | tr '|' '\n' | tr '\' '\n' | sed 's/CLNDBN[=]//g' | sort -nr  | head -100 > ${TOP_100}

awk '{$1=""}1' ${TOP_100} | awk '{$1=$1}1' | sort | uniq > ${TEMP_FILE}
mv ${TEMP_FILE} ${TOP_100}
