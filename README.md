# PSST
Polygenic SNP Search Tool

## Overview

Here is what this pipeline does: We are developing a software *pipeline* to **identify multiple SNPs that are associated with diseases**; including SNPs that modify the *penetrance* of other SNPs. We are looking to identify:
* Asserted pathogenic SNPs
* Genome-wide Association Studies (GWAS) identified SNPs
, crossed with database and datasets such as ClinVar, SRA, and GEO, and then build a pipeline for multiple genetic variants associated with diseases.


## Usage

### Step 1 -- Extracting phenotypes

#### This script extracts 100 disease phenotypes that have the most associated SNPs in ClinVar that are asserted to be pathogenic

##### This truncates extensions and details and takes root terms

###### It does not do any semantic clustering of like phenotypes

### Step 2 -- Generating a summary table

#### This script makes a table that summarizes the number of SNPs and pathogenic variants associated with the 100 disease phenotypes, as well as the number of SRA and GEO datasets available that are related to those phenotypes

### Step 3 -- Specifying SNPs and Datasets for selected phenotypes 

