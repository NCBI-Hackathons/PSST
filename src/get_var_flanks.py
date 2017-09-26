#!/usr/bin/env python
from Bio import Entrez
import sys
import getopt

def get_var_flanking_sequences(accessions,email):
    flanking_sequences = {}
    Entrez.email = email
    for var_id in accessions:
        var_id = var_id.rstrip()
        if var_id.startswith('rs'):
            var_id = var_id[len('rs'):] 
        var_id = var_id.rstrip()
        if len(var_id) > 0:
            handle = Entrez.esummary(db='snp',id=var_id,retmode='xml') 
            records = Entrez.parse(handle)
            for record in records:
                docsum = record['DOCSUM']
                docsum_tokens = docsum.split('|')
                flanking_seq = [token for token in docsum_tokens if 'SEQ=' in token][0].split('=')[1]
                flanking_sequences[var_id] = flanking_seq
            handle.close()
    return flanking_sequences

def write_flanking_sequences(flanking_sequences,output_path):
    with open(output_path,'w') as out_stream:
        for var_id in flanking_sequences:
            sequence = flanking_sequences[var_id]
            line = "%s=%s\n" % (var_id,sequence) 
            out_stream.write(line)

def unit_tests():
    accessions = ['187525243']
    print(get_var_flanking_sequences(accessions,'laseanl@sfu.ca'))

def main(argv):
    help_message = "Description: given a list of variant accessions, gets their flanking sequences and outputs\n"\
                 + "             them into the specified output file."
    usage_message = "Usage: [-i variant accessions file] [-e email address for Entrez servers] [-o output path]"
    options = "hti:e:o:"

    try:
        opts,args = getopt.getopt(argv[1:],options)
    except getopt.GetoptError:
        print("Error: unable to read command line arguments.")
        sys.exit(1)

    if len(argv) == 1:
        print(help_message)
        print(usage_message)
        sys.exit(0)

    accessions_file = None
    output_path = None
    email = None

    for opt, arg in opts:
        if opt == '-h':
            print(help_message)
            print(usage_message)
        elif opt == '-i':
            accessions_file = arg
        elif opt == '-e':
            email = arg
        elif opt == '-o':
            output_path = arg
        elif opt == '-t':
            unit_tests()
            sys.exit(0)

    opts_incomplete = False

    if accessions_file == None:
        print("Error: please provide the path to the accessions file.")
        opts_incomplete = True
    if output_path == None:
        print("Error please provide a path for the output file.")
        opts_incomplete = True
    if opts_incomplete:
        sys.exit(1)

    accessions = []

    with open(accessions_file,'r') as in_stream:
        for line in in_stream:
            accessions.append( line.rstrip() )     
    flanking_sequences = get_var_flanking_sequences(accessions,email)
    write_flanking_sequences(flanking_sequences,output_path)

if __name__ == "__main__":
    main(sys.argv)
