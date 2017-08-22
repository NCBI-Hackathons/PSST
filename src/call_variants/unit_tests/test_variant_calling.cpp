#include <string>
#include <vector>
#include <map>
#include "catch.hpp"
#include "../structs.hpp"
#include "../variant_calling.hpp"

TEST_CASE( "call_variants returns the correct variants", "[call_variants]" )
{
    std::vector<std::string> variant_accessions;
    std::map<std::string,VariantFrequency> map;
    VariantFrequency variant_freq;

    std::string accession = "het";
    variant_freq.matches = 4;
    variant_freq.mismatches = 6;
    variant_accessions.push_back(accession);
    map[accession] = variant_freq; 

    accession = "hom";
    variant_freq.matches = 10;
    variant_freq.mismatches = 0;
    variant_accessions.push_back(accession);
    map[accession] = variant_freq; 

    accession = "nil";
    variant_freq.matches = 0;
    variant_freq.mismatches = 10;
    variant_accessions.push_back(accession);
    map[accession] = variant_freq; 

    VariantFrequenciesMap variant_frequencies_map;
    variant_frequencies_map.variant_accessions = variant_accessions;
    variant_frequencies_map.map = map;

    CalledVariants called_variants = call_variants(variant_frequencies_map);

    REQUIRE( called_variants.homozygous.size() == 1 );
    REQUIRE( called_variants.heterozygous.size() == 1 );
}

TEST_CASE( "alignment_spans_variants returns the correct values", "[alignment_spans_variant]" )
{
    VarBoundary var_boundary;
    var_boundary.start = 1;
    var_boundary.stop = 2;
    Alignment alignment;
    // alignment does not span variant
    alignment.ref_start = 3;
    alignment.ref_stop = 10;
    REQUIRE( not alignment_spans_variant(alignment,var_boundary) );
    // alignment does span variant
    alignment.ref_start = 0;
    REQUIRE( alignment_spans_variant(alignment,var_boundary) );
}
