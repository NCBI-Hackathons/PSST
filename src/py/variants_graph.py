#!/usr/bin/env python
## Copyright: NCBI 2017
## Authors: Sean La, Emmanuel Marte
import getopt
import sys
import networkx as nx
from itertools import combinations
import matplotlib.pyplot as plt

def read_variant_tsv(tsv_path):
    '''
    Reads the TSV file outputted by the `call_variants` executable in the PSST repository.
    Outputs a dictionary where the keys are SRA accessions and the return value is another dictionary containing the
    homozygous and heterozygous variants. 
    '''
    sra_variants = {}
    with open(tsv_path,'r') as tsv:
        for line_no, line in enumerate(tsv):
                if line_no > 0:
                    tokens = line.split('\t')
                    sra_accession = tokens[0]
                    heterozygous = tokens[1].split(' ')
                    homozygous = tokens[2].split(' ') 
                    variants = {'heterozygous':heterozygous,'homozygous':homozygous}
                    sra_variants[sra_accession] = variants
    return sra_variants

def construct_graph(sra_variants):
    '''
    Given the variant dictionary outputted by `read_variant_tsv`, constructs a graph of the related variants
    ''' 
    variant_graph = nx.Graph()
    for sra_accession in sra_variants:
        variants = sra_variants[sra_accession]
        heterozygous = variants['heterozygous']
        homozygous = variants['homozygous']
        all_variants = heterozygous + homozygous
        pairs = list(combinations(all_variants,2))
        for pair in pairs:
            if not variant_graph.has_edge(*pair):
                variant_graph.add_edge(*pair,weight=1)
            else:
                u = pair[0]
                y = pair[1]
                variant_graph[u][v]['weight'] += 1
    return variant_graph 

if __name__ == '__main__':
    help_message =  \
    """
    Displays related variants as a graph structure where variants are vertices and SRA datasets are edges.
    Two vertices are connected by an edge if and only if they both appear in the same SRA dataset.
    The weight of an edge uv is the number of SRA datasets that u and v both appear in.
    """ 

    usage_message = "%s [-h help and usage] [-i input variant TSV file] [-o output path for graph]" % (sys.argv[0])

    options = "hi:o:"

    try:
        opts, args = getopt.getopt(sys.argv[1:],options)
    except getopt.GetoptError:
        print("Error: unable to read command line arguments.")
        sys.exit(1)

    if len(sys.argv) == 1:
        print(usage_message)
        sys.exit(0)

    tsv_path = None
    output_path = None

    for opt,arg in opts:
        if opt == '-h':
            print(help_message)
            print(usage_message)
            sys.exit(0)
        elif opt == '-i':
            tsv_path = arg
        elif opt == '-o':
            output_path = arg

    opts_incomplete = False

    if tsv_path == None:
        opts_incomplete = True
        print("Error: please provide the variant TSV file.")
    if output_path == None:
        opts_incomplete = True
        print("Error: please specify an output path.")
    if opts_incomplete:
        print(usage_message)
        sys.exit(1)

    sra_variants = read_variant_tsv(tsv_path)
    variant_graph = construct_graph(sra_variants)
    nx.draw(variant_graph)
    plt.show()
