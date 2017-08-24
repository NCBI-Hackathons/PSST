#!/usr/bin/env python
# Copyright: NCBI 2017
# Authors: Sean La
import getopt
import sys
from get_alleles import get_nth_allele

help_message = "Given a file containing lines of the form 'ACCESSION=W[X/Y]Z', this script creates a FASTA file\n" \
             + "where the header identifiers are 'ACCESSION' and the sequence is the major allele, i.e. 'WXZ'" 
usage_message = "[-h help and usage] [-i flanking sequence file] [-o output file]"

options = "hi:o:"

try:
	opts, args = getopt.getopt(sys.argv[1:],options)
except getopt.GetoptError:
	print("Error: unable to read command line arguments.")
	sys.exit(1)

if len(sys.argv) == 1:
	print(usage_message)
	sys.exit(0)

input_path = None
output_path = None

for opt, arg in opts:
	if opt == '-h':
		print(help_message)
		print(usage_message)
		sys.exit(0)
	elif opt == '-i':
		input_path = arg
	elif opt == '-o':
		output_path = arg

opts_incomplete = False

if input_path == None:
	print("Error: please provide the path to the input flank file.")
	opts_incomplete = True
if output_path == None:
	print("Error: please provide the path to the output FASTA file.")
	opts_incomplete = True

if opts_incomplete:
	print(usage_message)
	sys.exit(1)

with open(input_path,'r') as input_stream, open(output_path,'w') as output_stream:
	for line in input_stream:
		tokens = line.split('=')
		if len(tokens) == 2:
			accession = tokens[0]
			sequence = tokens[1]
			allele = get_nth_allele(sequence,2)
			output_stream.write( ">%s\n" % (accession) )
			output_stream.write( "%s" % (allele) )
