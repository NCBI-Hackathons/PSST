#!/usr/bin/env python
import sys
import argparse
from get_alleles import get_nth_allele

def find_var_info(sequences):
	'''
	Finds the length of the flanking sequences left and right of a variant.
	Also returns the length of the minor allele
	Input
	- sequences: dictionary of sequences
	Output
	- var_info: dictionary of flanking lengths for each sequence
	'''
	var_info = {}
	for seq_name in sequences:
		sequence = sequences[seq_name]
		minor_allele = get_nth_allele(sequence,2)

		left_bracket_index = sequence.find('[')
		right_bracket_index = sequence.find(']')
		start = left_bracket_index

		right_flank_length = len(sequence[right_bracket_index+1:])
		length = len(minor_allele)
		stop = length - right_flank_length
		var_info[seq_name] = (start,stop,length)
	return var_info

def unit_test():
	# Set up test
	seq_name = "test"
	sequence = "AAAA[GA/T]AAAAA"
	major_allele = "AAAAGAAAAAA"
	real_length = len(major_allele)
	real_left = 4
	right_right = 5

	sequences = { seq_name : sequence }
	var_info = find_var_info(sequences)
	flanks = var_info[seq_name]
	start = flanks[0]
	stop = flanks[1]
	length = flanks[2]
	assert(start == 4)
	assert(stop == 5) 
	assert(length == 11)
	print("Unit tests passed!")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(add_help=False,description=
	'''
	Author: Sean La. Given a list of sequences through STDIN in the form "SEQ_NAME=W[X/Y]Z" where W,X,Y, and Z\
	are nucleotide subsequences, prints a line in the form "SEQ_NAME START STOP LENGTH" where START and STOP\
        are the start and stop indices of the variant X and LENGTH is the length of the allele containing X. 
	''')
	parser.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,
				help='Show this help message and exit.')
	parser.add_argument('-t','--test',action='store_true',help=
	"""
	Perform unit tests for this script.
	""")
	parser.add_argument('-i','--input',metavar='INPUT',help=
	"""
	Path to the input file. If not set, then takes input through STDIN.
	""")
	parser.add_argument('-o','--output',metavar='OUTPUT',help=
	"""
	Path to the output file. If not set, outputs through STDOUT.
	""")
	args = parser.parse_args()

	if args.test:
		unit_test()
		sys.exit(0)
	
	if args.input:
		input_stream = open(args.input,'r')
	else:
		input_stream = sys.stdin
	
	sequences = {} 

	for line in input_stream:
		tokens = line.split("=")
		if len(tokens) == 2:
			name = tokens[0]
			sequence = tokens[1]
			sequences[name] = sequence

	input_stream.close()

	var_info = find_var_info(sequences)

	info_lines = []
	for seq_name in var_info:
		start = var_info[seq_name][0]
		stop = var_info[seq_name][1]
		length = var_info[seq_name][2]
		info_lines.append( "%s %d %d %d" % (seq_name, start, stop, length) )

	if args.output:
		with open(args.output,'w') as output_stream:
			for line in info_lines:
				output_stream.write( "%s\n" % (line) )
	else:
		for line in info_lines:
			print(line) 
