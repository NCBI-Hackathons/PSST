#!/usr/bin/env python
import sys
import argparse

# Global variables are depicted in all uppercase
MATCH = '=' # Base to represent a match in the alignment
GAP = '_' # Base to represent a gap in the alignment
QUERY = 0 # The first base of a BTOP base-pair represents the variant in the query
REF = 1 # The second base of a BTOP base-pair represents the variant in the reference

def find_delimited_btop(btop):
	'''
	Constructs a BTOP string so that there is a space between each operation
	Example: 4CGAT_10_4 -> 4 CG AT _10_ 4
	Inputs
	- (str) btop: the BTOP string
	Outputs
	- (str): The delimited BTOP string
	'''
	btop.replace("^","_") # We treat introns as gaps
	delimited_btop = ""
	num_base = 0
	in_gap = False
	prev_was_int = False
	for i in range( len(btop) ):
		current_char = btop[i]
		if current_char.isalpha() or current_char == '-':
			num_base += 1
			if num_base % 2 == 1:
				delimited_btop += " "
			prev_was_int = False
		elif current_char == "_":
			if in_gap:
				in_gap = False
			else:
				delimited_btop += " "
				in_gap = True
			prev_was_int = False	
		elif current_char.isdigit():
			if not prev_was_int and not in_gap:
				delimited_btop += " "
			prev_was_int = True	
		delimited_btop += current_char
	return delimited_btop

def delimited_btop_to_alignment(delimited_btop):
	'''
	Given a delimited BTOP string, outputs the reference sequence alignment 
	Inputs
	- (str) delimited_btop: the space delimited BTOP string
	Outputs
	- (str) ref: the reference sequence alignment
	'''
	# Construct an alignment out of the delimited BTOP string 
	btop_list = delimited_btop.split()
	ref = ""
	gap = False
	for btop in btop_list:
		if '_' in btop:
			gap = True
			btop = ''.join(btop.split('_'))
		else:
			gap = False	
		if btop.isdigit():
			num_matches = int(btop)
			if gap:
				base = GAP
			else:
				base = MATCH
			ref = ref + (base * num_matches)
		else:
			ref = ref + btop[REF]	
	return ref

def read_flank_info(path):
	'''
	Given the path to a file containing the flank information for a set of variants, returns a dict with that info
	Inputs
	- (str) path to the file
	Outputs
	- a dictionary where the keys are the variant IDs and the values are a list containing the flank information
	'''
	flank_info = {}
	with open(path,'r') as flanks:
		for line in flanks:
			tokens = line.split()
			# The first item is the ID of the SNP
			# The rest are the flank lengths and length of the major allele sequence
			flank_info[tokens[0]] = {'left':tokens[1],'right':tokens[2],'length':tokens[3]}
	return flank_info
			

def translate_var_boundary(left,ref):
	'''
	Translates the left boundary of the variant of interest in the reference sequence into the left boundary
	of the variant in the reference alignment
	e.g.
	Given ref = ATCG, left = 1 and ref_alignment = A-TCG, this script outputs alignment_left = 2
	Inputs
	- (int) left: the length of the left flank in the alignment
	- (str) ref: the reference alignment
	Outputs
	- (int) alignment_left: the length of the right most boundary of the left flank in the alignment
	'''
	# find the position of the ref base in the reference alignment
	alignment_left = 0 # This will be the position of the ref base in the reference alignment
	seq_left = 0 # Keeps track of the current position in the reference sequence 
	while seq_left < left - 1 and alignment_left - 1 < len(ref):
		current_base = ref[alignment_left]
		if current_base != "-":
			seq_left += 1
		alignment_left += 1
	return alignment_left + 1

def query_contains_ref_base(alignment,flank_info):
	'''
	Determines whether the query sequence as encoded by the BTOP string contains the ref base at pos.
	Inputs
	- alignments: contains information regarding the magic-blast alignment for the query
	- flank_info: lengths of the SNP flanks and length of the major allele
	Outputs
	- (bool): True if query contains ref bases surrounded by the flanks, False otherwise 
	'''
	# Extract the information in the alignment and flank_info objects
	btop = alignment['btop']
	ref_start = alignment['ref_start']
	ref_stop = alignment['ref_stop']
	left = flank_info['left']
	right = flank_info['right']	
	length = flank_info['length']
	# Determines whether the aligned subsequence of the reference even includes the variant
	if ref_start > left or ref_stop < right: 
		return False
	else:
		# Adjust the variant boundary indices to be compatible with the local alignment given by the BTOP
		# string
		left = left - ref_start
		right = right - ref_start
	# Convert the BTOP string into an alignment
	delimited_btop = find_delimited_btop(btop)
	# Get the reference sequence alignment
	ref = delimited_btop_to_alignment(delimited_btop)
	# Translate the variant boundaries in the sequence into the boundaries for the alignment 
	alignment_left = translate_var_boundary(left,ref)
	alignment_right = translate_var_boundary(right,ref)
	# Start and stop positions for the interval containing the variant in the alignment
	start = alignment_left - 1
	stop = len(ref) - alignment_right
	# If any base in the interval containing the variant is not a match, then the query does not contain the
	# variant
	for i in range(start,stop):
		ref_base = ref[i]
		if ref_base != MATCH:
			return False
	return True 

def unit_tests():
	# Set up test data
	orig_ref = "AAAAGTTTTTTTTTTAAAA"
	orig_query = "AAAACCAAAA"
	btop = "4C-CG_10_4"
	left = 16
	right = 17
	alignment = {'btop':btop,'ref_start':0,'ref_stop':len(orig_ref)}
	flank_info = {'left':left,'right':right,'length':len(orig_ref)}

	# Start unit tests
	delimited_btop = find_delimited_btop(btop)
	assert(delimited_btop == " 4 C- CG _10_ 4")

	ref = delimited_btop_to_alignment(delimited_btop)
	assert(ref == "====-G__________====")

	alignment_left = translate_var_boundary(left,ref)
	alignment_right = translate_var_boundary(right,ref)
	assert(alignment_left == left + 1)
	assert(alignment_right == right + 1)

	assert( query_contains_ref_base(alignment,flank_info) )

	alignment = {'btop':btop, 'ref_start':0, 'ref_stop':len(orig_ref)}
	flank_info = {'left':5, 'right':6, 'length':len(orig_ref)}
	assert( not query_contains_ref_base(alignment,flank_info) )

	print("All unit tests passed!")

if __name__ == '__main__':
	# Columns in Magic-BLAST tabulated output for the query name and BTOP string
	query_k = 0
	btop_k = 16
	# Set up the argument parser
	parser = argparse.ArgumentParser(add_help=False,description=
	'''
	Author: Sean La. Given Magic-BLAST tabulated output provided either through STDIN or as a file,\
	returns a list of queries that contain the same base as the reference at a specified location in the \
	reference. Indexing of sequences is 0-based. 
	''',
	epilog="Example: magic-blast.tab | queriesWithRefBase.py 21 >> queries_with_variant.txt")
	parser.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,
                    help='Show this help message and exit.')
	parser.add_argument('-i','--input_file',metavar='INPUT_FILE',type=str,help=
		"""
		Path to file containing Magic-BLAST tabulated output. If not set, will read through STDIN.	
		""")
	parser.add_argument('-f','--flank_info',metavar='FLANK_INFO',type=str,help=
		"""
		Path to the file containing information regarding the lengths of the variant flanking sequences
		""")
	parser.add_argument('-o','--output_path',metavar='OUTPUT_PATH',type=str,help=
		"""
		Output path for list of queries that contain the ref base. If not set, will output query names in 
		STDOUT.
		""")
	parser.add_argument('-t','--test',action="store_true",help=
		"""
		Perform unit tests for this script.
		""")
	parser.add_argument('left',type=int,help='Length of left flanking sequeuence in reference sequence.')
	parser.add_argument('right',type=int,help='Length of right flanking sequence in reference sequence.')
	args = parser.parse_args()
	# Perform unit tests and exit if specified
	if args.test:
		unit_tests()
		sys.exit(0)
	# If the user provided an input file, read from it. Otherwise, read from stdin
	if args.input_file:
		input_stream = open(args.input_file,'r')
	else:
		input_stream = sys.stdin
	# Find those queries that contain the specified reference base
	queries = []
	for line in input_stream:
		if len(line) > 0 and line[0] != '#':
			tokens = line.split()
			btop = tokens[btop_k]
			if query_contains_ref_base(btop,args.left,args.right):
				query = tokens[query_k]
				queries.append(query)
	input_stream.close()
	# If the user provided an output file path, write the queries there. Otherwise, just print to stdout.
	if args.output_path:
		with open(args.output_path,'w') as output_stream:
			for query in queries:
				output_stream.write( "%s\n" % (query) )
	else:
		for query in queries:
			print(query)		
