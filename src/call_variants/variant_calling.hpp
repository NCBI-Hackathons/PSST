#ifndef VARIANT_CALLING_H
#define VARIANT_CALLING_H

#include "structs.hpp"

CalledVariants call_variants(VariantFrequenciesMap var_freq_map);

bool alignment_spans_variant(Alignment alignment, VarBoundary var_boundary);

#endif /* VARIANT_CALLING_H */
