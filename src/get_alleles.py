def get_major_allele(seq):
        '''
        Given a variant sequence in the form 'W[X/Y]Z' where W,X,Y,Z are nucleotide sequences, returns the major\
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
        # Construct the major allele
        major_allele = seq[:left_bracket_index] + major_var + seq[right_bracket_index + 1:]
        return major_allele
