#include <vector>
#include <string>
#include "catch.hpp"
#include "../data_manipulation.hpp"
#include "../structs.hpp"

TEST_CASE( "split returns a vector whose size is one more than the number of spaces in the input string", "[split]")
{
    std::string test_string = "aaaa bbbb cccc";
    int num_tabs = 2;
    std::vector< std::string > tokens = split(test_string,' ');
    REQUIRE( tokens.size() == num_tabs + 1 );
}

TEST_CASE("split returns the correct vector given a string","[split]")
{
    std::string test_string = "aaaa bbbb cccc";
    std::vector<std::string> split_tokens = split(test_string,' ');
    std::vector<std::string> actual_tokens;
    actual_tokens.push_back("aaaa");
    actual_tokens.push_back("bbbb");
    actual_tokens.push_back("cccc");
    for (int i = 0; i < actual_tokens.size(); i++) {REQUIRE(split_tokens.at(i) == actual_tokens.at(i));}
}
