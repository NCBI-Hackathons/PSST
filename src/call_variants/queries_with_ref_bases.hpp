#ifndef QUERY_WITH_REF_BASES
#define QUERY_WITH_REF_BASES

#include "structs.hpp"

bool is_number(const std::string &str);

std::vector<std::string> get_btop_list(std::string btop);
/* Given a BTOP alignment string (outputted by Magic-BLAST), returns a vector where each item is an individual BTOP
 * operation. */

std::string find_reference_alignment(std::vector<std::string> btop_list);
/* Given a BTOP list, returns a pseudoalignment of the reference where matches are not explicitly expressed */ 

int64_t get_alignment_var_boundary(std::string ref, int64_t boundary);
/* Given the boundary of a variant in the reference sequence, returns the boundary of the variant in the reference
 * alignment. */

bool query_contains_ref_bases(Alignment alignment, VarBoundary var_boundary);
/* Determines whether the given query contains the same sequence of bases at a designated area as the reference in
 * the alignment. */

#endif /* QUERY_WITH_REF_BASES */
