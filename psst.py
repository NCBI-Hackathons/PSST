#!/usr/bin/env python
## Built-in packages
import getopt
import sys
import os
import subprocess
import shlex
## Custom packages
from src import get_var_flanks
from src import find_var_info
from src import var_flanks_to_fasta
from src import call_variants
from src import get_alleles

## Command line arguments
MAJOR_VERSION = 2
MINOR_VERSOIN = 1

parser = argparse.ArgumentParser(description='''
    Long Read Correction Stats (LRCstats) Version %d.%d Copyright (C) 2017 Sean La
    This program comes with ABSOLUTELY NO WARRANTY
    This is free software, and you are welcome to redistribute it
    under certain conditions
    ''' % (MAJOR_VERSION, MINOR_VERSION))

parser.add_argument('-s','--sra', metavar='SRA_FILE', type=str, help="Path to the file containing SRA accessions")
parser.add_argument('-n','--snp', metavar='SNP_FILE', type=str, help="Path to the file containing SNP accessions")
parser.add_argument('-f','--fastq', metavar='FASTQ', type=str, help="Path to the NGS FASTQ file")
parser.add_argument('-d','--dir', metavar='WORKING_DIR', type=str, help="Directory to store temporary files")
parser.add_argument('-e','--email', metavar='EMAIL', type=str, help="Email address for Entrez")
parser.add_argument('-t','--threads', metavar='THREADS', type=int, help="Number of threads for each Magic-BLAST call")
parser.add_argument('-c','--child_procs', metavar='CHILD_PROCS', type=int, help="Maximum number of child processes")

args = parser.parse.args()

## Make sure the arguments are complete
opts_incomplete = False 

if args.sra and args.fastq:
    opts_incomplete = True
    print("Error: please provide only one of SRA accessions or FASTQ file, not both.")
elif (not args.sra) and (not args.fastq):
    opts_incomplete = True
    print("Error: please provide either an SRA accessions file or a FASTQ file.")
if not args.snp:
    opts_incomplete = True
    print("Error: please provide a SNP accessions file.")
if not args.dir:
    opts_incomplete = True
    print("Error: please specify a working directory.")
if not args.email:
    opts_incomplete = True
    print("Error: please provide an email address for Entrez.")
if not args.threads:
    opts_incomplete = True
    print("Error: please specify the number of threads to be used with each Magic-BLAST run.")
if not args.child_procs:
    opts_incomplete = True
    print("Error: please specify the maximum number of child processes.")

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

# Get the directory that this script exists in, i.e. the path to the PSST directory on the machine
src = os.path.realpath(__file__) + "/src" 

## Create BLAST database out of the variant FASTA file
print("Creating a BLAST database out of the SNP flanking sequences...")
makeblastdb_command = "%s/makeblastdb.sh %s %s" % (src, fasta_path, args.dir)
subprocess.call( shlex.split(makeblastdb_command) )

## Align the SRA dataset(s) or FASTQ file onto the SNP flanking sequences using Magic-BLAST
if args.sra:
    blast_command = "%s/magicblast_sra.sh %s snp_flanks %s %d %d" \
                    % (src,args.sra,args.dir,args.threads,args.child_procs)
else:
    blast_command = "%s/magicblast_fastq.sh %s snp_flanks %s %d" % (src,args.fastq,args.dir,args.threads) 
subprocess.call( shlex.split(blast_command) )

## Call variants
print("Calling variants...")
accession_map = {str( snp_acc.index(accession) ): accession for accession in snp_acc}
mbo_directory = args.dir + "/mbo" 
paths = call_variants.get_mbo_paths(mbo_directory)
map_paths_and_partition = {'map':accession_map,'paths':paths,'partition':paths.keys()}
sra_alignments = call_variants.get_sra_alignments(map_paths_and_partition)
alignments_and_info = {'alignments':sra_alignments,'keys':sra_alignments.keys(),'info':var_info}
called_variants = call_variants.call_sra_variants(alignments_and_info)
output_path = src + "/results.tsv"
create_tsv(called_variants,output_path)
