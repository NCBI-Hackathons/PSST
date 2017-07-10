#!/usr/bin/env python
from Bio import Entrez
import sys
import argparse
import subprocess

def get_major_allele(seq):
	'''
	Given a known variant in the form 'SEQ=W[X/Y]Z' where W,X,Y,Z are nucleotide sequences, returns the major\
	allele, i.e. WXZ
	Input
	- seq: the variant sequence
	Output
	- the major allele as a string
	'''
	left_bracket_index = seq.find('[')
	right_bracket_index = seq.find(']')	
	# extract the two variants
	variants = seq[ left_bracket_index + 1 : right_bracket_index ]
	var_tokens = variants.split('/')
	major_var = var_tokens[0]
	minor_var = var_tokens[1]
	# Construct the major allele
	major_allele = seq[4:left_bracket_index] + major_var + seq[right_bracket_index + 1:]
	return major_allele

def write_fasta(variants,temp):
	'''
	Given known variants, writes a FASTA file in temp
	Inputs
	- (dict) variants: the RNAseq variants
	- (str) temp: directory to contain the FASTA file
	Returns
	- (dict) paths: paths to the FASTA files
	'''
	paths = {}
	for varId in variants:
		fasta_path = "%s/%s.fasta" % (temp,varId)
		with open(fasta_path,'w') as fasta:
			fasta.write( ">%s\n" % (varId) )
			fasta.write( "%s\n" % (variants[varId]) ) 
		paths[varId] = fasta_path
	return paths

def make_magicblastdb(fasta_paths,temp):
	'''
	Given a dictionary of RNAseq variant FASTA files, creates magicblast databases
	Inputs
	- (dict) fasta_paths: dictionary containing paths to fasta files
	- (str) temp: directory to contain magicblast db files	
	'''
	for varId in fasta_paths:
		fasta_path = fasta_paths[varId]
		command = "makeblastdb -dbtype nucl -in %s -out %s/%s" % (fasta_path,temp,varId)
		subprocess.Popen(command)

# Arguments parser
parser = argparse.ArgumentParser(add_help=False,description=
		'''
		Author: Sean La. Given the SEQ ID of a known variant in RNA-seq data, does stuff with it.
		''')
parser.add_argument('email',type=str,help='Email to make use of Entrez')
parser.add_argument('id',type=str,help='ID of the variant')
parser.add_argument('temp',type=str,help='Directory to hold temporary files')
parser.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit.')
args = parser.parse_args()

# Always tell NCBI who you are
Entrez.email = args.email
# Find the data from the SNP database
Entrez.tool = "Polygenic Variant Checker"
search_handle = Entrez.esearch(db="snp",term=args.id)
search_record = Entrez.read(search_handle)
# For now, just take the first item in the IdList. We might use all of them in the future.
summary_handle = Entrez.esummary(db="snp",id=search_record["IdList"][0],retmax=search_record["RetMax"])
# For some reason Entrez.read returns a dictionary within a list, so we take out the dict out of the list
summary_record = Entrez.read(summary_handle)[0]
# Get the DOCSUM info
docsum = summary_record["DOCSUM"]
docsum_tokens = docsum.split('|')
# Extract the sequence from DOCSUM
seq = [token for token in docsum_tokens if "SEQ" in token]
# We assume there is only one sequence for the variant 
seq = seq[0]
# Take the major allele
major_allele = get_major_allele(seq)
# Combine all this information into one neat package
variants = { args.id : major_allele }
# Create FASTA files for the variant
paths = write_fasta(variants, args.temp)
# Construct a magicblast database out of the variants
make_magicblastdb(paths,args.temp)
