#!/bin/bash
# Copyright: NCBI 2017
# Author: Chipo Mashayomombe and Sean La
 
if [ "$#" -ne 2 ]; then
	echo "Description: Generates a table that summarizes the number of SNPs and pathogenic variants associated"
	echo "             with the top 100 disease phenotypes as well as the number of related SRA and GEO datasets."
	BASENAME=`basename "$0"`
	echo "Usage: ${BASENAME} [Top 100 phenotypes file] [Output path]"
	exit 0
fi

TOP100=$1
OUTPUT=$2

## Create the columns of the table

DISEASE=""
PATHOGENIC=""  
SNPs=""
SRA=""
GEO=""

for PHENOTYPE in $(cat ${TOP100})
do
	DISEASE="${DISEASE} \t $PHENOTYPE"
	# number of pathogenic variants associated with the phenotype
	PATHOGENIC="${PATHOGENIC} \t $(esearch -query $PHENOTYPE -db clinvar|efilter -query "pathogenic" | xtract -pattern ENTREZ_DIRECT -element Count)"
	# number of snps associated with the phenotype 
	SNP="${SNP} \t $(esearch -query $PHENOTYPE -db pubmed | elink -target snp | esummary | xtract -pattern DocumentSummary -element SNP_ID | wc -l)"
	# number of SRA datasets associated with the disease phenotype
	SRA="${SRA} \t $(esearch -query $PHENOTYPE -db sra |xtract -pattern ENTREZ_DIRECT -element Count)"
	# number of GEO datasets associated with the disease phenotype
	GEO="${GEO} \t $(esearch -query $PHENOTYPE -db gds |xtract -pattern ENTREZ_DIRECT -element Count)"
done


## Create the summary table file
echo "" > ${OUTPUT}

for LINE in DISEASE PATHOGENIC SNP SRA GEO
do
	echo -e ${LINE} >> ${OUTPUT}
done

echo "Summary table construction complete."
