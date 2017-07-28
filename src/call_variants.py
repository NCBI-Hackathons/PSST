#!/usr/bin/env python
# Built-in python packages
from __future__ import division # This modifies Python 2.7 so that any expression of the form int/int returns a float 
import getopt
import sys
import os
# Project-specific packages
from queries_with_ref_bases import query_contains_ref_bases

def get_accession_map(fasta_path):
	'''
	Suppose in the FASTA file used as reference for makeblastdb, there are n sequences.
	Magic-BLAST labels these sequences as the order in which they appear when outputting in tabulated format.
	For example, suppose that the SNP rs0001 appears second from the top in the FASTA file.
	Then Magic-BLAST assigns the label '1' to rs0001 whenever it appears in an alignment.
	This function returns a dictionary that serves as a map from integers (in string datatype) to accessions
	e.g. accession_map['1'] == rs0001 in our above example.
	Inputs
	- (str) fasta_path: path to the FASTA file used as reference for makeblastdb
	Outputs
	- (dict) accession_map: the map from integers to accessions
	'''
        accession_map = {}
        with open(fasta_path,'r') as fasta:
                id_number = 0
                for line in fasta:
                        if line[0] == ">":
                                accession = line[1:].rstrip()
                                accession_map[str(id_number)] = accession
                                id_number += 1
        return accession_map


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
			# We only want the accession number, not the extension as well
			accession = os.path.basename(file).split('.')[0]
			path = os.path.join(directory,file)
			paths[accession] = path
	return paths

def get_sra_alignments(paths,accession_map):
	'''
	Given a list of paths as described in the function get_mbo_paths, retrieves the BTOP string for each
	alignment.
	Inputs
	- paths: a list of pairs where the first entry of the pair is the accession and the second is the path
	- (dict) accession_map: the map between integers and accessions
	Outputs
	- a dictionary where keys are SRA accessions and the values are alignment dictionaries
	'''
	sra_alignments = {}
	for accession in paths:
		path = paths[accession]
		alignments = []	
		with open(path,'r') as mbo:
			for line in mbo:
				tokens = line.split()
				# Skip the line if it is commented, the number of fields isn't equal to 25 or
				# the query read was not aligned
				if line[0] != "#" and len(tokens) == 25 and tokens[1] != "-":
					var_acc = accession_map[ tokens[1] ]
					ref_start = tokens[8]
					ref_stop = tokens[9]
					if int(ref_start) > int(ref_stop):
						temp = ref_start
						ref_start = ref_stop
						ref_stop = temp
					btop = tokens[16]
					alignment = { 'var_acc': var_acc, 'ref_start': ref_start,\
						      'ref_stop': ref_stop, 'btop': btop }
					alignments.append( alignment )
		sra_alignments[accession] = alignments
	return sra_alignments

def get_var_info(path):
	'''
	Retrieves the flanking sequence lengths for the SNP sequences			
	Inputs
	- path: path to the file that contains the flanking sequence lengths
	Outputs
	- var_info: a dictionary where the keys are SNP accessions and the values are tuples where the first entry\
                    is the start position of the variant, the second is the stop position of the variant and \
                    the third is the length of the variant sequence
	'''
	var_info = {}
	with open(path,'r') as input_file:
		for line in input_file:
			tokens = line.split()
			if len(tokens) == 4:
				accession = tokens[0]
				start = int(tokens[1])
				stop = int(tokens[2])
				length = int(tokens[3])
				var_info[accession] = {'start':start,'stop':stop,'length':length}
	return var_info

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

def get_sra_variants(sra_alignments,var_info):
	'''
	For all SRA accession, determines which variants exist in the SRA dataset	
	Inputs
	- sra_alignments: dict where the keys are SRA accessions and the values are lists of alignment dicts
	- var_info: dict where the keys are variant accessions and the values are information concerning the variants
	Outputs
	- variants: dict where the keys are the SRA accessions and the values are lists which contain the accessions
		    of variants which exist in the SRA datasets
	'''
	variants = {}
	for sra_acc in sra_alignments:
		alignments = sra_alignments[sra_acc]
		var_freq = {}
		for alignment in alignments:
			var_acc = alignment['var_acc']
			# Get the flank information
			info = var_info[var_acc]
			if var_acc not in var_freq:
				var_freq[var_acc] = {'true':0,'false':0}
			# Determine whether the variant exists in the particular SRA dataset
			var_called = query_contains_ref_bases(alignment,info)
			if var_called:
				var_freq[var_acc]['true'] += 1	
			else:
				var_freq[var_acc]['false'] += 1	
		sra_variants = call_variants(var_freq) 
		variants[sra_acc] = sra_variants	
	return variants

def create_tsv(variants,output_path):
	'''
	Creates a TSV file containing the set of variants each SRA dataset contains.
	Inputs
	- variants: dict where the keys are the SRA accessions and the values are lists which contain the accessions
		    of variants which exist in the SRA datasets
	- output_path: path to where to construct the output file
	'''
	with open(output_path,'w') as tsv:
		header = "SRA	Variants\n"
		tsv.write(header)
		for sra_acc in variants:
			line = "%s" % (sra_acc)		
			sra_variants = variants[sra_acc]
			for var_acc in sra_variants:
				line = "%s	%s" % (line, var_acc)
			tsv.write(line)
			tsv.write('\n')

if __name__ == "__main__":
	help_message = "Description: Given a directory with Magic-BLAST output files where each output file\n" \
                     + "             contains the alignment between an SRA dataset and known variants in a human\n" \
                     + "             genome, this script determines which variants each SRA dataset contains\n" \
                     + "             using a heuristic."
	usage_message = "Usage: %s [-h (help and usage)]\n[-m <directory containing .mbo files>]\n" % (sys.argv[0]) \
                      + "              [-v <path to variant info file>]\n[-f <path to the reference FASTA file]\n"\
                      + "              [-o <output path for TSV file>]"
	options = "hm:v:f:o:"

	try:
		opts,args = getopt.getopt(sys.argv[1:],options)
	except getopt.GetoptError:
		print("Error: unable to read command line arguments.")
		sys.exit(1)

	if len(sys.argv) == 1:
		print(help_message)
		print(usage_message)
		sys.exit()

	mbo_directory = None
	var_info_path = None
	output_path = None
	fasta_path = None
	
	for opt, arg in opts:
		if opt == '-h':
			print(help_message)
			print(usage_message)
			sys.exit(0)
		elif opt == '-m':
			mbo_directory = arg
		elif opt == '-v':
			var_info_path = arg
		elif opt == '-o':
			output_path = arg
		elif opt == '-f':
			fasta_path = arg
		elif opt == '-t':
			unit_tests()
			sys.exit(0)

	opts_incomplete = False

	if mbo_directory == None:
		print("Error: please provide the directory containing your Magic-BLAST output files.")
		opts_incomplete = True
	if var_info_path == None:
		print("Error: please provide the path to the file containing flanking sequence information.")
		opts_incomplete = True
	if output_path == None:
		print("Error: please provide an output path for the TSV file.")
		opts_incomplete = True
	if fasta_path == None:
		print("Error: please provide the path to the FASTA file used as reference for makeblastdb")
		opts_incomplete = True
	if opts_incomplete:
		print(usage_message)
		sys.exit(1)

	paths = get_mbo_paths(mbo_directory)
	accession_map = get_accession_map(fasta_path)
	sra_alignments = get_sra_alignments(paths,accession_map)
	var_info = get_var_info(var_info_path)
	variants = get_sra_variants(sra_alignments,var_info)
	create_tsv(variants,output_path)
