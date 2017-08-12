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
            flank_info[tokens[0]] = {'start':int(tokens[1]),'stop':int(tokens[2]),'length':int(tokens[3])}
    return flank_info
            

def translate_var_boundary(boundary,ref):
    '''
    Translates the start boundary of the variant of interest in the reference sequence into the start boundary
    of the variant in the reference alignment
    e.g.
    Given ref = ATCG, start = 1 and ref_alignment = A-TCG, this script outputs alignment_start = 2
    Inputs
    - (int) boundary: the index of the boundary of the variant
    - (str) ref: the reference alignment
    Outputs
    - (int) alignment_boundary: the length of the stop most boundary of the start flank in the alignment
    '''
    # find the position of the ref base in the reference alignment
    alignment_boundary = 0 # This will be the position of the ref base in the reference alignment
    seq_index = 0 # Keeps track of the current position in the reference sequence 
    while seq_index < boundary and alignment_boundary < len(ref):
        current_base = ref[alignment_boundary]
        if current_base != "-":
            seq_index += 1
        alignment_boundary += 1
    return alignment_boundary

def query_contains_ref_bases(alignment,flank_info):
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
    start = flank_info['start']
    stop = flank_info['stop']    
    length = flank_info['length']
    # Determines whether the aligned subsequence of the reference even includes the variant
    if ref_start > start or ref_stop < stop: 
        return None
    else:
        # Adjust the variant boundary indices to be compatible with the local alignment given by the BTOP
        # string
        start = start - ref_start
        stop = stop - ref_start
    # Convert the BTOP string into an alignment
    delimited_btop = find_delimited_btop(btop)
    # Get the reference sequence alignment
    ref = delimited_btop_to_alignment(delimited_btop)
    # Translate the variant boundaries in the sequence into the boundaries for the alignment 
    alignment_start = translate_var_boundary(start,ref)
    alignment_stop = translate_var_boundary(stop,ref)
    # Start and stop positions for the interval containing the variant in the alignment
    start = alignment_start - 1
    stop = len(ref) - alignment_stop
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
    start = 16
    stop = 17
    alignment = {'btop':btop,'ref_start':0,'ref_stop':len(orig_ref)}
    flank_info = {'start':start,'stop':stop,'length':len(orig_ref)}

    # Start unit tests
    delimited_btop = find_delimited_btop(btop)
    assert(delimited_btop == " 4 C- CG _10_ 4")

    ref = delimited_btop_to_alignment(delimited_btop)
    assert(ref == "====-G__________====")

    alignment_start = translate_var_boundary(start,ref)
    assert(alignment_start == start + 1)
    alignment_stop = translate_var_boundary(stop,ref)
    assert(alignment_stop == stop + 1)

    assert( query_contains_ref_bases(alignment,flank_info) )

    alignment = {'btop':btop, 'ref_start':0, 'ref_stop':len(orig_ref)}
    flank_info = {'start':5, 'stop':6, 'length':len(orig_ref)}
    assert( not query_contains_ref_bases(alignment,flank_info) )

    print("All unit tests passed!")

if __name__ == '__main__':
    # Set up the argument parser
    parser = argparse.ArgumentParser(add_help=False,description=
    '''
    Author: Sean La. Given Magic-BLAST tabulated output information, this script determines which query\
        sequences contain the same bases as the reference sequences.
    ''')
    parser.add_argument('-h','--help',action='help',default=argparse.SUPPRESS,
                    help='Show this help message and exit.')
    parser.add_argument('-t','--test',action="store_true",help=
        """
        Perform unit tests for this script.
        """)
    args = parser.parse_args()
    # Perform unit tests and exit if specified
    if args.test:
        unit_tests()
        sys.exit(0)
