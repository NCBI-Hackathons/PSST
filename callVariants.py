#!/usr/bin/env python
import argparse
import os
from __future__ import division # This modifies Python 2.7 so that any expression of the form int/int returns a float 

# Global variables with accessing elements within tuples
ACC = 0
PATH = 1

def get_mbo_paths(directory):
	'''
	Given a directory, retrieves the paths of all files within the directory whose extension is .mbo
	Inputs
	- (str) directory: directory to be scoured for mbo files
	Outputs
	- paths: a dictionary where the keys are accessions and the values are paths
	'''
	paths = {}
	for file in os.listdir(directory):
		if file.endswith(".mbo"):
			accession = os.path.basename(file)	
			path = os.path.join(directory,file)
			paths[accession] = path
	return paths

def get_btops(paths):
	'''
	Given a list of paths as described in the function get_mbo_paths, retrieves the BTOP string for each
	alignment.
	Inputs
	- paths: a list of pairs where the first entry of the pair is the accession and the second is the path
	Outputs
	- a dictionary where keys are SRA accessions and the values are alignment dictionaries
	'''
	btops = {}
	for accession in paths:
		path = pair[accession]
		alignments = []	
		with open(path,'r') as mbo:
			for line in mbo:
				tokens = line.split()
				snp_acc = tokens[1]
				ref_start = tokens[8]
				ref_stop = tokens[9]
				btop = tokens[16]
				alignment = { 'snp_acc': tokens[1], 'ref_start': tokens[8],\
					      'ref_stop': tokens[9], 'btop': tokens[16] }
				alignments.append( alignment )
		btops[accession] = alignments
	return btops

def get_flanks(path):
	'''
	Retrieves the flanking sequence lengths for the SNP sequences			
	Inputs
	- path: path to the file that contains the flanking sequence lengths
	Outputs
	- flanks: a dictionary where the keys are SNP accessions and the values are pairs where the first item
                  is the length of the left flank and the second is the length of the right flank
	'''
	flanks = {}
	with open(path,'r') as file:
		for line in file:
			tokens = line.split()
			accession = tokens[0]
			left = tokens[1]
			right = tokens[2]
			flanks[accession] = (left,right)	
	return flanks

def call_snps(snp_freq):
	'''
	Determines which snps exist in a given SRA dataset given the number of reads that do and do not contain
	the snp
	Inputs
	- snp_freq: a dict where the keys are SNP accessions and the values are dicts which contain the frequency of
		    reads that do and reads that do not contain the SNP
	Outputs
	- snps: a list which contains the SNP accessions of those SNPs which exist in the SRA dataset
	'''
	snps = []
	for snp_acc in snp_freq:
		frequencies = snp_freq
		true = frequencies['true']
		false = frequencies['false']
		percentage = true / (true + false)
		if percentage > 0.4: # For now, we use this simple heuristic.
			snps.append(snp_acc)
	return snps

def get_sra_snps(btops,flanks)
	'''
	For all SRA accession, determines which snps exist in the SRA dataset	
	Inputs
	- btops: dict where the keys are SRA accessions and the values are lists of alignment dicts
	- flanks: dict where the keys are SNP accessions and the values are the lengths of the flanks
	Outputs
	- snps: dict where the keys are the SRA accessions and the value is a list which contains the accessions
		    of SNPs which exist in the SRA dataset
	'''
	snps = {}
	for sra_acc in btops:
		snp_freq = {}
		alignments = btops[sra_acc]
		for alignment in alignments:
			snp_acc = alignments['snp_acc']
			if snp_acc not in snp_freq:
				snp_freq[snp_acc] = { 'true': 0, 'false': 0 }
			snp_called = query_contains_ref_base(alignment)
			if snp_called:
				snp_freq[snp_acc]['true'] += 1	
			else:
				snp_freq[snp_acc]['false'] += 1	
		sra_snps = call_snps(snp_freq) 
		snps[sra_acc] = sra_snps	
	return snps
			

parser = argparse.ArgumentParser(description='')
parser.add_argument('directory', help=
	"""
	Input directory to the Magic-BLAST output files in tabulated format.
	""")
parser.add_argument('flanks', help=
	"""
	Input file for lengths of flanks for each SNP sequence
	""")
args = parser.parse_args()

paths = get_mbo_paths(args.directory)
btops = get_btops(paths)
flanks = get_flanks(args.flanks)
snps = get_sra_snps(btops,flanks)
