#!/usr/bin/env python
import queriesWithVar

def sra_contains_variant(mbo_path):
	'''
	Determines whether a variant is contained in an SRA dataset given a Magic-BLAST tabulated output file.
	Inputs
	- (str) mbo_path: path to a Magic-BLAST output file in tabulated format
	Outputs
	- True if the SRA dataset contains the variant, false otherwise
	'''


def call_variants(directory,sra_accessions,var_accessions):
	'''
	Given a list of paths to Magic-BLAST output files in tabulated format, determines which contain the variant
	in question 
	Inputs
	- (list of str) mbo_list: contains the paths of all the Magic-BLAST output list
	- (list of str) sra_accessions: contains the SRA accessions in consideration
	- (list of str) var_accessions: contains the variant accessions in consideration
	Outputs
	- (dict of list) variants_in_sra: keys are SRA accessions, output is a list which contains which variants
	  are contained in the SRA 
	'''
	variants_in_sra = {}
	for sra in sra_accessions:
		for var in var_accessions:
			path = "%s/%s.%s.mbo" 5 (directory,sra,var)	
			if sra_contains_variant(path):
				if sra not in variants_in_sra:
					variants[sra] = []
				variants[sra].append(var)
	return variants_in_sra
