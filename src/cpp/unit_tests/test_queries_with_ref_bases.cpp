#include <iostream>
#include <string>
#include <vector>
#include "catch.hpp"
#include "../queries_with_ref_bases.hpp"


TEST_CASE( "get_btop_list operates correctly", "[get_btop_list]" )
{
    std::string btop = "4C-CG_10_4";
    std::vector<std::string> btop_list = get_btop_list(btop);
    std::vector<std::string> expected_btop_list;
    expected_btop_list.push_back("4");
    expected_btop_list.push_back("C-");
    expected_btop_list.push_back("CG");
    expected_btop_list.push_back("_10_");
    expected_btop_list.push_back("4");
    REQUIRE( btop_list.size() == expected_btop_list.size() );
    for (int i = 0; i < expected_btop_list.size(); i++) { REQUIRE( btop_list.at(i) == expected_btop_list.at(i) ); }
}

TEST_CASE( "find_reference_alignment operates correctly", "[find_reference_alignment]" )
{
    std::string btop = "4C-CG_10_4";
    std::vector<std::string> btop_list = get_btop_list(btop);
    std::string reference_alignment = find_reference_alignment(btop_list);
    std::string actual_alignment = "====-G__________====";
    REQUIRE( reference_alignment.length() == actual_alignment.length() );
    REQUIRE( reference_alignment == actual_alignment );
}

TEST_CASE( "get_alignment_var_boundary works correctly", "[get_alignment_var_boundary]" )
{
    int64_t start = 4;
    int64_t stop = 5;
    std::string alignment = "====-G__________====";
    int64_t translated_start = get_alignment_var_boundary(alignment, start);
    int64_t translated_stop = get_alignment_var_boundary(alignment, stop);
    int64_t expected_start = 5;
    int64_t expected_stop = 6;
    REQUIRE( translated_start == expected_start );
    REQUIRE( translated_stop == expected_stop );

    start = 16;
    stop = 17;
    translated_start = get_alignment_var_boundary(alignment, start);
    translated_stop = get_alignment_var_boundary(alignment, stop);
    expected_start = 17;
    expected_stop = 18;
    REQUIRE( translated_start == expected_start );
    REQUIRE( translated_stop == expected_stop );
}

TEST_CASE( "query_contains_ref_bases works correctly", "[query_contains_ref_bases]" )
{
    Alignment alignment;
    alignment.variant_acc = "test";
    alignment.btop = "4C-CG_10_4";
    alignment.ref_start = 0;
    alignment.ref_stop = 10;
    VarBoundary var_boundary;
    var_boundary.start = 3;
    var_boundary.stop = 3;
    REQUIRE( query_contains_ref_bases(alignment,var_boundary) );
    var_boundary.start = 10;
    var_boundary.stop = 10;
    REQUIRE( not query_contains_ref_bases(alignment,var_boundary) );
}
