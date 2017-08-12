#!/usr/bin/env python
# Built-in python packages
from __future__ import division # This modifies Python 2.7 so that any expression of the form int/int returns a float 
import getopt
import sys
import os
from itertools import combinations
from multiprocessing.dummy import Pool
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

def get_sra_alignments(map_paths_and_partition):
    '''
    Given a list of paths as described in the function get_mbo_paths, retrieves the BTOP string for each
    alignment.
    Inputs
    - map_paths_and_partition: a dict which contains the following:
        - partition: the list of paths to .mbo files to read
        - paths: a list of pairs where the first entry of the pair is the accession and the second is the path
        - (dict) accession_map: the map between integers and accessions
    Outputs
    - a dictionary where keys are SRA accessions and the values are alignment dictionaries
    '''
    accession_map = map_paths_and_partition['map']
    paths = map_paths_and_partition['paths']
    partition = map_paths_and_partition['partition']
    sra_alignments = {}
    for accession in partition:
        path = paths[accession]
        alignments = []    
        with open(path,'r') as mbo:
            for line in mbo:
                tokens = line.split()
                # Skip the line if it is commented, the number of fields isn't equal to 25 or
                # the query read was not aligned
                if line[0] != "#" and len(tokens) == 25 and tokens[1] != "-":
                    var_acc = accession_map[ tokens[1] ]
                    ref_start = int(tokens[8])
                    ref_stop = int(tokens[9])
                    if ref_start > ref_stop:
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
    - variants: dict where the keys are SRA accessions and the value is another dict that contains the homozgyous and 
                heterozygous variants in separate lists 
    '''
    variants = {'heterozygous':[],'homozygous':[]}
    for var_acc in var_freq:
        frequencies = var_freq[var_acc]
        true = frequencies['true']
        false = frequencies['false']
        try:
            percentage = true/(true+false)
            if percentage > 0.8: # For now, we use this simple heuristic.
                variants['homozygous'].append(var_acc)
            elif percentage > 0.3:
                variants['heterozygous'].append(var_acc)
        except ZeroDivisionError: # We ignore division errors because they correspond to no mapped reads
            pass
    return variants

def call_sra_variants(alignments_and_info):
    '''
    For all SRA accession, determines which variants exist in the SRA dataset    
    Inputs
    - alignments_and_info: a dict which contains
        - sra_alignments: dict where the keys are SRA accessions and the values are lists of alignment dicts
        - var_info: dict where the keys are variant accessions and the values are information concerning the variants
        - keys: list which contains the keys of the SRA accessions to analyze
    Outputs
    - variants: dict where the keys are SRA accessions and the value is another dict that contains the homozgyous and 
                heterozygous variants in separate lists 
    '''
    sra_alignments = alignments_and_info['alignments']
    var_info = alignments_and_info['alignments']
    keys = alignments_and_info['keys']
    variants = {}
    for sra_acc in keys:
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
            if var_called == True:
                var_freq[var_acc]['true'] += 1    
            elif var_called == False:
                var_freq[var_acc]['false'] += 1    
        sra_variants = call_variants(var_freq) 
        variants[sra_acc] = sra_variants    
    return variants

def create_tsv(variants,output_path):
    '''
    Creates a TSV file containing the set of variants each SRA dataset contains.
    Inputs
    - variants: dict where the keys are SRA accessions and the value is another dict that contains the homozgyous and 
                heterozygous variants in separate lists 
    - output_path: path to where to construct the output file
    '''
    with open(output_path,'w') as tsv:
        header = "SRA\tHeterozygous SNPs\tHomozygous SNPs\n"
        tsv.write(header)
        for sra_acc in variants:
            line = "%s" % (sra_acc)        
            sra_variants = variants[sra_acc]
            line = line + "\t"
            for var_acc in sra_variants['heterozygous']:
                line = line + var_acc + ","
            line = line + "\t"
            for var_acc in sra_variants['homozygous']:
                line = line + var_acc + ","
            tsv.write( line.rstrip() )
            tsv.write('\n')

def create_variant_matrix(variants):
    '''
    Returns an adjacency matrix that represents the graph constructed out of the variants such that:
    1. Every vertex represents a variant and every variant is represented by a vertex
    2. An edge (line) connects two vertices if and only if there exists an SRA dataset that contains both of the
       corresponding variants
    The matrix is represented by a dictionary of dictionaries. To save memory, we do not store 0 entries or variants
    without any incident edges. 
    Inputs
    - variants: dict where the keys are SRA accessions and the value is another dict that contains the homozgyous and 
                heterozygous variants in separate lists 
    Outputs
    - matrix: a dict which, for any two keys (variants) variant_1 and variant_2, satisfies the following -
              1. type(matrix[variant_1]) is DictType
              2. type(matrix[variant_1][variant_2]) is IntType and matrix[variant_1][variant_2] >= 1
              3. matrix[variant_1][variant_2] == matrix[variant_2][variant_1] 
              all of which holds if and only if variant_1 and variant_2 exist in the dictionaries as keys
    '''
    matrix = {}
    for sra_acc in variants:
        sra_variants = variants[sra_acc]
        all_variants = []
        if 'homozygous' in sra_variants:
            all_variants += sra_variants['homozygous']
        if 'heterozygous' in sra_variants:
            all_variants += sra_variants['heterozygous']
        # Get all of the unique 2-combinations of the variants 
        two_combinations = list( combinations(all_variants,2) )
        for pair in two_combinations:
            variant_1 = pair[0]
            variant_2 = pair[1] 
            if variant_1 not in matrix:
                matrix[variant_1] = {}
            if variant_2 not in matrix:
                matrix[variant_2] = {}
            if variant_2 not in matrix[variant_1]: 
                matrix[variant_1][variant_2] = 0
            if variant_1 not in matrix[variant_2]:
                matrix[variant_2][variant_1] = 0
            matrix[variant_1][variant_2] += 1
            matrix[variant_2][variant_1] = matrix[variant_1][variant_2]
    return matrix

def partition(lst,n):
    '''
    Partitions a list into n lists
    Inputs
    - (list) lst
    Outputs
    - partitioned_lists: a list of lists 
    '''
    division = len(lst)/float(n)
    return [ lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n) ]

def combine_list_of_dicts(list_of_dicts):
    '''
    Given a list of dicts, returns a single dict such that (key,value) pairs from each original dict exists in the
    new single dict
    Inputs
    - list_of_dicts: a list of dicts
    Outputs
    - combined_dict 
    '''
    combined_dict = list_of_dicts[0]
    for i in range(1,len(list_of_dicts)):
        combined_dict.update( list_of_dicts[i] )
    return combined_dict

def unit_tests():
    variants = {}
    variants['sra_1'] = {'homozygous':['a','b'],'heterozygous':['c','e']}
    variants['sra_2'] = {'homozygous':['a','c','d']}
    variants['sra_3'] = {'heterozygous':['b','d','e']}
    matrix = create_variant_matrix(variants)
    print(matrix)
    for variant_1 in matrix:
        for variant_2 in matrix[variant_1]:
            left_hand_side = matrix[variant_1][variant_2]
            right_hand_side = matrix[variant_2][variant_1]
            assert( left_hand_side >= 1 )
            assert( right_hand_side >= 1 )
            assert( left_hand_side == right_hand_side )
    print("All unit tests passed!")

if __name__ == "__main__":
    help_message = "Description: Given a directory with Magic-BLAST output files where each output file\n" \
                     + "             contains the alignment between an SRA dataset and known variants in a human\n" \
                     + "             genome, this script determines which variants each SRA dataset contains\n" \
                     + "             using a heuristic."
    usage_message = "Usage: %s\n[-h (help and usage)]\n[-m <directory containing .mbo files>]\n" % (sys.argv[0]) \
                      + "[-v <path to variant info file>]\n[-f <path to the reference FASTA file>]\n"\
                      + "[-o <output path for TSV file>]\n[-p <num of threads>]\n[-t <unit tests>]"
    options = "htm:v:f:o:p:"

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
    threads = 1
    
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
        elif opt == '-p':
            threads = arg
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

    var_info = get_var_info(var_info_path)
    accession_map = get_accession_map(fasta_path)
    paths = get_mbo_paths(mbo_directory)

    # Retrieve the alignments concurrently
    get_alignments_threads = min(threads,len(paths.keys()))
    paths_partitions = partition( paths.keys(), get_alignments_threads )
    map_paths_and_partitions = [{'map':accession_map,'paths':paths,'partition':path_partition} \
                                for path_partition in paths_partitions]
    pool = Pool(processes=get_alignments_threads)
    sra_alignments_pool = pool.map(get_sra_alignments,map_paths_and_partitions)
    pool.close()
    pool.join()
    sra_alignments = combine_list_of_dicts(sra_alignments_pool) 

    # Call variants concurrently
    sra_keys = sra_alignments.keys()
    variant_call_threads = min( threads, len(sra_keys) )
    keys_partitions = partition(sra_keys, variant_call_threads)
    # alignments, info and key partitions
    alignments_and_info_part = [{'alignments':sra_alignments,'keys':keys,'info':var_info} for keys in keys_partitions]
    pool = Pool(processes=variant_call_threads)
    variants_pool = pool.map(call_sra_variants,alignments_and_info_part)
    pool.close()
    pool.join()
    called_variants = combine_list_of_dicts(variants_pool)

    create_tsv(called_variants,output_path)
    matrix = create_variant_matrix(called_variants)
