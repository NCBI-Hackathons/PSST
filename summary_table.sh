#!/bin/bash
# Copyright: NCBI 2017
# Author: Chipo Mashayomombe
 
# create variables to store the data for the disease phenotype as well as the
# number of 
#   -pathogenic variants associated with the phenotype (PATHOGENIC)
#   -snps associated with the phenotype (SNPs)
#   -SRA datasets associated with the disease phenotype (SRA)
#   -GEO datasets associated with the disease phenotype (GEO)

DISEASE=""
PATHOGENIC=""  
SNPs=""
SRA=""
GEO=""
 
# Assign the top 100 disease phenotypes to a variable
top100=$(cat top100_disease_phenotypes)
 
for PHENOTYPE in $top100
do
    DISEASE="${DISEASE} $PHENOTYPE"
    PATHOGENIC="${PATHOGENIC} $(esearch -query $PHENOTYPE -db clinvar|efilter -query "pathogenic" | xtract -pattern ENTREZ_DIRECT -element Count)"
    SNPs="${SNPs} $(esearch -query $PHENOTYPE -db pubmed | elink -target snp | esummary | xtract -pattern DocumentSummary -element SNP_ID | wc -l)"
    SRA="${SRA} $(esearch -query $PHENOTYPE -db sra |xtract -pattern ENTREZ_DIRECT -element Count)"
    GEO="${GEO} $(esearch -query $PHENOTYPE -db gds |xtract -pattern ENTREZ_DIRECT -element Count)"
done


