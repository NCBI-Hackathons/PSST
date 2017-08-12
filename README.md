# PSST
Polygenic SNP Search Tool Version 1.1

## Graphical Overview

![Workflow](/media/Polygenic_SNP_Search_Tool.png?raw=true "Workflow.png")

## Overview:

The Polygenic SNP Search Tool is an open-source pipeline that **identifies multiple SNPs that are associated with diseases**; including SNPs that modify the *penetrance* of other SNPs. This pipeline identifies:
* Asserted pathogenic SNPs
* Genome-wide Association Studies (GWAS) identified SNPs
, crossed with database and datasets such as ClinVar, SRA, and GEO, and then constructs a report describing multiple genetic variants associated with diseases.


## Usage:

The main script `psst.sh` accepts as input SRA (accession prefix `SRR`) and SNP (accession prefix `rs`) accessions where, in the files containing the accessions, each line corresponds to a unique accessions.
This script will then output a TSV file describing which SNPs are contained in the SRA datasets.

The `psst.sh` subpipeline is as follows:

1. Extracts flanking sequences for the SNP accessions and creates a FASTA file containing these flanking sequences. 

2. Uses `makeblastdb` to generate a BLAST database for the SNP flanking sequences.

3. Runs Magic-BLAST on each phenotype-associated SRA dataset and the SNP flanking sequence BLAST database.

4. From the Magic-BLAST alignments, determines which SNPs are contained in the SRA datasets using a statistical heuristic.

See the file `breast-ovarian_cancer.tsv` for an example output file.

## Disease Clustering:

Clustering disease types through the ClinVar database in various categories such as assorted metabolic diseases and breast cancer to see the relationship among human variations and phenotypes. 

1. Diseases were manually found using ClinVar dataset. Those that were not a match were eliminated. 
