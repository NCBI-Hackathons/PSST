# PSST
Polygenic SNP Search Tool

## Graphical Overview

![Workflow](/media/Polygenic_SNP_Search_Tool.png?raw=true "Workflow.png")

## Overview

Here is what this pipeline does: We are developing a software *pipeline* to **identify multiple SNPs that are associated with diseases**; including SNPs that modify the *penetrance* of other SNPs. We are looking to identify:
* Asserted pathogenic SNPs
* Genome-wide Association Studies (GWAS) identified SNPs
, crossed with database and datasets such as ClinVar, SRA, and GEO, and then build a pipeline for multiple genetic variants associated with diseases.


## Usage:

Step 1:

Generating a summary table:

* The script makes a table that summarizes the number of SNPs and pathogenic variants associated with the 100 disease phenotypes, as well as the number of SRA and GEO datasets available that are related to those phenotypes

Step 2:

The pipeline finds SNP IDs for disease phenotypes from a list of 100 disease phenotypes that have the SNPs associated in ClinVar with the disease and are asserted to be pathogenic.

Step 3:

Extracts flanking sequences for these SNP IDs, and also creates seperate fasta files for each disease containing both the SNP IDs and the respective sequences.

Step 4:

Uses makeblastdb to generate a database for each phenotype.

Step 5:

Performs alignment with magicblast using SRA accessions associated with disease phenotype.



