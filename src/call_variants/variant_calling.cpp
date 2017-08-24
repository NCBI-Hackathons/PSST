#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "structs.hpp"
#include "variant_calling.hpp"
#include "queries_with_ref_bases.hpp"

bool alignment_spans_variant(Alignment alignment, VarBoundary var_boundary)
/* Determines whether the aligned segment of the variant flanking sequence contains the variant */
{
    int64_t ref_start = alignment.ref_start;
    int64_t ref_stop = alignment.ref_stop;
    int64_t var_start = var_boundary.start;
    int64_t var_stop = var_boundary.stop;
    return (ref_start <= var_start and var_stop <= ref_stop);
}

CalledVariants call_variants_in_sra(VariantFrequenciesMap &var_freq_map)
/* Determine the set of variants that exist in a particular SRA dataset */
{   // Retrieve the data from the VariantFrequencies object
    std::vector<std::string> accessions = var_freq_map.variant_accessions;
    std::map<std::string,VariantFrequency> map = var_freq_map.map;
    // Initialize the output
    CalledVariants called_variants;
    std::vector<std::string> homozygous;
    std::vector<std::string> heterozygous;
    // Iterate through variant_frequencies
    for (size_t index = 0; index < accessions.size(); index++) {
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

std::vector<CalledVariants> call_variants(const std::vector<SraAlignments> &sra_alignments_list,
                                          std::map<std::string,VarBoundary> &var_boundary_map)
/* Determines which variants exist in multiple SRA datasets */
{
    std::vector<CalledVariants> called_variants_list;
    // Call variants in each of the SRA datasets
    for (size_t index = 0; index < sra_alignments_list.size(); index++) {
        SraAlignments sra_alignments = sra_alignments_list.at(index);
        VariantFrequenciesMap var_freq_map;
        // Analyze the alignments of the SRA dataset onto the variant flanking sequences
        for (size_t alignments_index = 0; alignments_index < sra_alignments.alignments.size(); alignments_index++) {
            Alignment alignment = sra_alignments.alignments.at(alignments_index);
            std::string variant_acc = alignment.variant_acc;
            VarBoundary var_boundary = var_boundary_map[variant_acc];
            // Make sure the alignment spans the variant first
            if ( alignment_spans_variant(alignment,var_boundary) ) {
                // Make sure the variant exists in var_freq_map
                if (var_freq_map.map.find(variant_acc) == var_freq_map.map.end()){
                    var_freq_map.variant_accessions.push_back(variant_acc);
                    var_freq_map.map[variant_acc].matches = 0;
                    var_freq_map.map[variant_acc].mismatches = 0;
                }
                //Determine whether there is a match between the read and the variant flanking sequence at the variant
                bool matched = query_contains_ref_bases(alignment,var_boundary);
                // Update the matches and mismatches counts
                if (matched) {var_freq_map.map[variant_acc].matches++;}
                else {var_freq_map.map[variant_acc].mismatches++;}
            }
        }
        // From all the variant match and mismatch counts, determine which variants exist in the SRA dataset
        CalledVariants called_variants = call_variants_in_sra(var_freq_map);
        called_variants.accession = sra_alignments.accession;
        called_variants_list.push_back(called_variants);
    }
    return called_variants_list;
}
