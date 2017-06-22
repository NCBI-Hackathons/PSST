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

def write_fasta(variants,output):
	'''
	Given known variants, writes a FASTA file in temp
	Inputs
	- (dict) variants: the RNAseq variants
	- (str) temp: directory to contain the FASTA file
	Returns
	- (dict) paths: paths to the FASTA files
	'''
	with open(output,'w') as fasta:
	# Add the flanking sequence for each SNP to the FASTA file
		for varId in variants:
			fasta.write( ">%s\n" % (varId) )
			fasta.write( "%s\n" % (variants[varId]) ) 

# Arguments parser
parser = argparse.ArgumentParser(add_help=False,description=
		'''
		Given a phenotype (e.g. "breast cancer"), returns the flanking sequences of all SNPs associated \
		with the input phenotype.
		''')
parser.add_argument('email',type=str,help='Email to make use of Entrez')
parser.add_argument('phenotype',type=str,help='Phenotype')
parser.add_argument('output',type=str,help='Output')
parser.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,help='Show this help message and exit.')
args = parser.parse_args()

# DEPRECRATED: do not use the code below for production.
'''
# Always tell NCBI who you are
Entrez.email = args.email
# Find the data from the SNP database
Entrez.tool = "Polygenic Variant Checker"
search_handle = Entrez.esearch(db="snp",term=args.output)
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
# Create FASTA file for the flanking sequences of the SNPs
paths = write_fasta(variants, args.temp)
'''
