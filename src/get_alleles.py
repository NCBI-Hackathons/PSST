def get_nth_allele(seq,n):
        '''
        Given a variant sequence in the form 'X[Y_1/../Y_t]Z' where X,Y_1,..,Y_t,Z are nucleotide sequences, returns
        the nth allele, i.e. X(Y_n)Z. If the nth allele does not exist, returns the last allele.
        Input
        - (str) seq: the variant sequence
	- (int) n: the desired allele
        Output
        - the nth allele as a string, or the last allele if the sequence has less than n alleles.
        '''
	# Find the boundaries of where the variants occur
        left_bracket_index = seq.find('[')
        right_bracket_index = seq.find(']')
        # extract the variants
        variants = seq[ left_bracket_index + 1 : right_bracket_index ]
        var_tokens = variants.split('/')
	num_alleles = len(var_tokens)
	# Choose either the nth or the last allele
	index = min([num_alleles,n]) - 1
        variant = var_tokens[index]
        # Construct the allele
        allele = seq[:left_bracket_index] + variant + seq[right_bracket_index + 1:]
        return allele
