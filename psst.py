#!/usr/bin/env python
## Built-in packages
import argparse
import sys
import os
## Custom packages
from src import get_var_flanks
from src import find_var_info
from src import var_flanks_to_fasta
from src import call_variants
from src import get_alleles

## Command line arguments
MAJOR_VERSION = 1
MINOR_VERSOIN = 0

parser = argparse.ArgumentParser(description=
    '''
    Polygenic SNP Search Tool (PSST) Version %d.%d Copyright (C) 2017 NCBI
    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software, and you are welcome to redistribute it
    under certain conditions.
    ''' % (MAJOR_VERSION, MINOR_VERSION) )
parser.add_argument('sra', metavar='SRA_FILE', type=str, help="Path to the file containing SRA accessions")
parser.add_argument('snp', metavar='SNP_FILE', type=str, help="Path to the file containing SNP accessions")
parser.add_argument('dir', metavar='WORKING_DIR', type=str, help="Directory to store temporary files")
parser.add_argument('email', metavar='EMAIL', type=str, help="Email address for Entrez")
parser.add_argument('threads', metavar='THREADS', type=int, help="Number of threads for each Magic-BLAST call")
parser.add_argument('child_procs', metavar='CHILD_PROCS', type=int, help="Maximum number of child processes")

args = parser.parse.args()

## Retrieve SRA and SNP accessions
print("Reading SRA and SNP accessions files...")
sra_acc = []
snp_acc = []
with open(args.sra,'r') as sra_file:
    for line in sra_file:
        sra_acc.append( line.rstrip() )
with open(args.snp,'r') as snp_file:
    for line in snp_file:
        snp_acc.append( line.rstrip() )

## Create the working directory if it doesn't exist yet
print("Creating working directory...")
os.makedirs(args.dir)

## Find the variant flanking sequences
print("Finding SNP flanking sequences...")
flanking_sequences = get_var_flanks.get_var_flanking_sequences(snp_acc,args.email)

## Get the variant flank length information
print("Getting SNP flank length information")
var_info = find_var_info.find_var_info(flanking_sequences)

## Construct a variant FASTA file
print("Creating SNP flank FASTA file...") 
fasta_path = "%s/snp_flanks.fasta" % (args.dir)
with open(fasta_path,'w') as fasta:
    for var_id in flanking_sequences:
        flanking_sequence = flanking_sequences[var_id]
        allele = get_alleles.get_nth_allele(flanking_sequence,2).rstrip() # Take the minor allele
        fasta.write( ">%s\n" % (var_id) )
        fasta.write( "%s\n" % (allele) )

## Create BLAST database out of the variant FASTA file

