#ifndef QUERY_WITH_REF_BASES
#define QUERY_WITH_REF_BASES

std::string find_delimited_btop(std::string btop);

std::string find_reference_alignment(std::string delimited_btop);

int64_t translate_var_boundary(std::string ref, int64_t boundary);

bool query_contains_ref_bases(AlignmentInfo alignment_info, FlankInfo flank_info);

#endif /* QUERY_WITH_REF_BASES */
