#ifndef VARIANT_CALLING_H
#define VARIANT_CALLING_H

#include <vector>
#include "structs.hpp"

bool alignment_spans_variant(Alignment alignment, VarBoundary var_boundary);

CalledVariants call_variants_in_sra(VariantFrequenciesMap &var_freq_map);

std::vector<CalledVariants> call_variants(const std::vector<SraAlignments> &sra_alignments_list,
                                          std::map<std::string,VarBoundary> &var_boundary_map);

#endif /* VARIANT_CALLING_H */
