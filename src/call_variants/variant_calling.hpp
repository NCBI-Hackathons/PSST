#ifndef VARIANT_CALLING_H
#define VARIANT_CALLING_H

#include <vector>
#include "structs.hpp"

bool alignment_spans_variant(Alignment alignment, VarBoundary var_boundary);
/* Determines whether the variant described in var_boundary is spanned by the alignment. */

CalledVariants call_variants_in_sra(VariantFrequenciesMap &var_freq_map);
/* Determines which variants exist in a particular SRA dataset by its variant frequencies. */

std::vector<CalledVariants> call_variants(const std::vector<SraAlignments> &sra_alignments_list,
                                          std::map<std::string,VarBoundary> &var_boundary_map);
/* Determines which variants exist in multiple SRA datasets. */

#endif /* VARIANT_CALLING_H */
