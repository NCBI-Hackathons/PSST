#!/usr/bin/env python
# Built-in python packages
from __future__ import division # This modifies Python 2.7 so that any expression of the form int/int returns a float 
import getopt
import sys
import os
# Project-specific packages
from query_with_ref_base import query_contains_ref_bases

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

def get_sra_alignments(paths):
	'''
	Given a list of paths as described in the function get_mbo_paths, retrieves the BTOP string for each
	alignment.
	Inputs
	- paths: a list of pairs where the first entry of the pair is the accession and the second is the path
	Outputs
	- a dictionary where keys are SRA accessions and the values are alignment dictionaries
	'''
	sra_alignments = {}
	for accession in paths:
		path = pair[accession]
		alignments = []	
		with open(path,'r') as mbo:
			for line in mbo:
				tokens = line.split()
				var_acc = tokens[1]
				ref_start = tokens[8]
				ref_stop = tokens[9]
				btop = tokens[16]
				alignment = { 'var_acc': tokens[1], 'ref_start': tokens[8],\
					      'ref_stop': tokens[9], 'btop': tokens[16] }
				alignments.append( alignment )
		sra_alignments[accession] = alignments
	return sra_alignments

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
	with open(path,'r') as input_file:
		for line in input_file:
			tokens = line.split()
			accession = tokens[0]
			left = tokens[1]
			right = tokens[2]
			length = tokens[3]
			flanks[accession] = {'left':left,'right':right,'length':length}
	return flanks

def call_variants(var_freq):
	'''
	Determines which variants exist in a given SRA dataset given the number of reads that do and do not contain
	the var
	Inputs
	- var_freq: a dict where the keys are SNP accessions and the values are dicts which contain the frequency of
		    reads that do and reads that do not contain the SNP
	Outputs
	- variants: a list which contains the SNP accessions of those SNPs which exist in the SRA dataset
	'''
	variants = []
	for var_acc in var_freq:
		frequencies = var_freq[var_acc]
		true = frequencies['true']
		false = frequencies['false']
		percentage = true/(true+false)
		if percentage > 0.4: # For now, we use this simple heuristic.
			variants.append(var_acc)
	return variants

def get_sra_variants(sra_alignments,flanks):
	'''
	For all SRA accession, determines which variants exist in the SRA dataset	
	Inputs
	- sra_alignments: dict where the keys are SRA accessions and the values are lists of alignment dicts
	- flanks: dict where the keys are SNP accessions and the values are the lengths of the flanks
	Outputs
	- variants: dict where the keys are the SRA accessions and the value is a list which contains the accessions
		    of SNPs which exist in the SRA dataset
	'''
	variants = {}
	for sra_acc in sra_alignments:
		alignments = sra_alignments[sra_acc]
		var_freq = {}
		for alignment in alignments:
			var_acc = alignment['var_acc']
			var_flanks = flanks[var_acc]
			if var_acc not in var_freq:
				var_freq[var_acc] = {'true':0,'false':0}
			# Determine whether the variant exists in the particular SRA dataset
			var_called = query_contains_ref_bases(alignment,var_flanks)
			if var_called:
				var_freq[var_acc]['true'] += 1	
			else:
				var_freq[var_acc]['false'] += 1	
		sra_variants = call_variants(var_freq) 
		variants[var_acc] = sra_variants	
	return variants

if __name__ == "__main__":
	help_message = ""
	usage_message = "%s [-h (help and usage)] [-m <directory containing .mbo files>] [-f <path to flanks info file>]" % (sys.argv[0])
	options = "hm:f:"
	try:
		opts,args = getopt.getopt(sys.argv[1:],options)
	except getopt.GetoptError:
		print("Error: unable to read command line arguments.")
		sys.exit(1)

	if len(sys.argv) == 1:
		print(usageMessage)
		sys.exit()

	mbo_directory = None
	flanks_path = None
	
	for opt, arg in opts:
		if opt == '-h':
			print(help_message)
			print(usage_message)
			sys.exit(0)
		elif opt == '-m':
			mbo_directory = arg
		elif opt == '-f':
			flanks_path = arg
		elif opt == '-t':
			unit_tests()
			sys.exit(0)

	optsIncomplete = False

	if mbo_directory == None:
		print("Error: please provide the directory containing your Magic-BLAST output files.")
		optsIncomplete = True
	if flanks_path == None:
		print("Error: please provide the path to the file containing flanking sequence information.")
		optsIncomplete = True
	if optsIncomplete:
		print(usage_message)
		sys.exit(1)

	paths = get_mbo_paths(mbo_directory)
	sra_alignments = get_sra_alignments(paths)
	flanks = get_flanks(flanks_path)
	variants = get_sra_variants(sra_alignments,flanks)
	print(variants)
