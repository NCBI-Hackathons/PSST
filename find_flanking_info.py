#!/usr/bin/env python
import sys
import argparse

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

def find_flank_info(sequences):
	'''
	Finds the length of the flanking sequences left and right of a variant.
	Assumes there is only one variant in the sequence.
	Input
	- sequences: dictionary of sequences
	Output
	- flanking_info: dictionary of flanking lengths for each sequence
	'''
	flanking_info = {}
	for seq_name in sequences:
		sequence = sequences[seq_name]
		left_bracket_index = sequence.find('[')
		right_bracket_index = sequence.find(']')
		left_flank_length = left_bracket_index
		right_flank_length = len(sequence[right_bracket_index+1:])
		length = len( get_major_allele(sequence) )
		flanking_info[seq_name] = (left_flank_length,right_flank_length,length)
	return flanking_info

def unit_test():
	sequence = "AAAA[G/T]AAAA"
	seq_name = "test"
	sequences = { seq_name : sequence }
	flanking_info = find_flank_info(sequences)
	flanks = flanking_info[seq_name]
	left_flank = flanks[0]
	right_flank = flanks[1]
	assert(left_flank == 4)
	assert(right_flank == 4) 
	print("Unit tests passed!")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(add_help=False,description=
	'''
	Author: Sean La. Given a list of sequences through STDIN in the form "seq_name=W[X/Y]Z" where W,X,Y, and Z\
	are nucleotide subsequences, returns the length of flanks W and Z through STDOUT in the form "seq_name left\
	right", where left and right are the lengths of the flanking sequences. Assumes each sequence has only\
	one variant.
	''')
	parser.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,
				help='Show this help message and exit.')
	parser.add_argument('-t','--test',action='store_true',help=
	"""
	Perform unit tests for this script.
	""")
	args = parser.parse_args()

	if args.test:
		unit_test()
		sys.exit(0)
	
	sequences = {} 
	input_stream = sys.stdin
	
	for line in input_stream:
		tokens = line.split("=")
		name = tokens[0]
		sequence = tokens[1]
		sequences[name] = sequence

	input_stream.close()

	flanking_info = find_flank_info(sequences)

	for seq_name in flanking_info:
		left_flank_length = flanking_info[seq_name][0]
		right_flank_length = flanking_info[seq_name][1]
		length = flanking_info[seq_name][2]
		print( "%s %d %d %d" % (seq_name, left_flank_length, right_flank_length, length) )
