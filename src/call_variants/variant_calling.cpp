#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "structs.hpp"
#include "variant_calling.hpp"

bool alignment_spans_variant(Alignment alignment, VarBoundary var_boundary)
/* Determines whether the aligned segment of the variant flanking sequence contains the variant */
{
    int64_t ref_start = alignment.ref_start;
    int64_t ref_stop = alignment.ref_stop;
    int64_t var_start = var_boundary.start;
    int64_t var_stop = var_boundary.stop;
    return (ref_start <= var_start and var_stop <= ref_stop);
}

CalledVariants call_variants(VariantFrequenciesMap var_freq_map)
/* Determine the set of variants that exist in a particular SRA dataset */
{   // Retrieve the data from the VariantFrequencies object
    std::vector<std::string> accessions = var_freq_map.variant_accessions;
    std::map<std::string,VariantFrequency> map = var_freq_map.map;
    // Initialize the output
    CalledVariants called_variants;
    std::vector<std::string> homozygous;
    std::vector<std::string> heterozygous;
    // Iterate through variant_frequencies
    for (int64_t index = 0; index < accessions.size(); index++) {
        std::string accession = accessions.at(index);
        VariantFrequency mapped_reads = map[accession];
        int64_t matches = mapped_reads.matches;
        int64_t mismatches = mapped_reads.mismatches;
        int64_t total_reads = matches + mismatches;
        if (total_reads > 0) {
            float percent_match = (float)matches/(float)total_reads;
            // We use this simple heuristic for now
            if (percent_match > 0.8) {homozygous.push_back(accession);}
            else if (percent_match > 0.3) {heterozygous.push_back(accession);}
        }
    }
    called_variants.homozygous = homozygous;
    called_variants.heterozygous = heterozygous;
    return called_variants;
}

