# PSST
Polygenic SNP Search Tool Version 2.0

## Graphical Overview

![Workflow](/media/Polygenic_SNP_Search_Tool4.png?raw=true "Workflow.png")

## Overview:

The Polygenic SNP Search Tool is an open-source pipeline that **identifies multiple SNPs that are associated with diseases**; including SNPs that modify the *penetrance* of other SNPs. This pipeline identifies:
* Asserted pathogenic SNPs
* Genome-wide Association Studies (GWAS) identified SNPs
, crossed with database and datasets such as ClinVar, SRA, and GEO, and then constructs a report describing multiple genetic variants associated with diseases.


## Usage:

The main script `psst.sh` accepts as input a text file where each line corresponds to a unique SNP rs-accessions and either another text file containing unique SRA accessions or a FASTQ file.
This script will then output a TSV file describing which SNPs are contained in the SRA datasets.

The PSST pipeline is as follows:

1. Extracts flanking sequences for the SNP accessions and creates a FASTA file containing these flanking sequences. 

2. Creates a BLAST database out of the SNP flanking sequences.

3. Runs Magic-BLAST on each phenotype-associated SRA dataset and the SNP flanking sequence BLAST database.

4. From the Magic-BLAST alignments, determines which SNPs are contained in the SRA datasets using a statistical heuristic.

See the file `breast-ovarian_cancer.tsv` for an example output file.

## Disease Clustering:

Grouping different disease types through the ClinVar database in various categories such as assorted metabolic diseases and breast cancer to see the relationship among human variations and phenotypes. 

1. Diseases were manually found exploring through the ClinVar dataset. 

2. Performed an online search to crosscheck whether the diseases that came up were metabolic or cancer related. 

3. Those that were not a match were eliminated while the correct diseases were moved into another file. 

## Future Additions
* Add a Bayesian inference variant calling rule for small number of NGS datasets. Our current heuristic runs fast on a large number of datasets, but for small number of datasets, a bayesian inference rule would be better and we wouldn't lose much in terms of time usage.
