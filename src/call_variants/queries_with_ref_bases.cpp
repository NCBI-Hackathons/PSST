#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>
#include <vector>
#include "structs.hpp"
#include "data_manipulation.hpp" // for split
#include "queries_with_ref_bases.hpp"

// Global variables
char GAP_g = '_';
char MATCH_g = '=';
int REF_g = 1;

bool is_number(const std::string &str)
{
    std::string::const_iterator it = str.begin();
    while (it != str.end() and std::isdigit(*it)) {++it;}
    return (not str.empty() and it == str.end());
}

std::vector<std::string> get_btop_list(std::string btop) 
{   // We treat introns as gaps
    std::replace( btop.begin(), btop.end(), '^', '_' );
    std::string temp_btop = "";
    int num_base = 0;
    bool in_gap = false;
    bool prev_was_int = false;
    for (size_t i = 0; i < btop.length(); i++) {
        char current_char = btop[i];
        if (std::isalpha(current_char) or current_char == '-') {
            num_base ++;
            if (num_base % 2 == 1) { temp_btop = temp_btop + " "; }
            prev_was_int = false;
        } else if (current_char == '_') {
            if (in_gap) { in_gap = false; }
            else { 
                temp_btop = temp_btop + " ";
                in_gap = true;
            }
            prev_was_int = false;
        } else if (std::isdigit(current_char)) {
            if (not prev_was_int and not in_gap) { temp_btop = temp_btop + " "; }
            prev_was_int = true;
        }
        temp_btop = temp_btop + current_char;
    }
    std::vector<std::string> btop_list = split(temp_btop);
    return btop_list;
}

std::string find_reference_alignment(std::vector<std::string> btop_list)
{ 
    std::string ref_alignment = "";
    std::string base; 
    bool gap = false;
    for (size_t i = 0; i < btop_list.size(); i++) {
        std::string op = btop_list.at(i);
        if (op.find("_") != std::string::npos) {
            gap = true;
            op.erase(std::remove(op.begin(),op.end(),'_'),op.end());
        } else {
            gap = false;
        }
        if ( is_number(op) ) {
            int num_matches = std::stoi(op);
            if (gap) { base = GAP_g; } else { base = MATCH_g; }
            for (int u = 0; u < num_matches; u++) { ref_alignment = ref_alignment + base; }
        } else { ref_alignment = ref_alignment + op[REF_g]; }
    }
    return ref_alignment;
}

int64_t get_alignment_var_boundary(std::string ref_alignment, int64_t boundary) 
{ 
    int64_t alignment_index = -1;
    int64_t seq_index = -1;
    int64_t ref_length = static_cast<int>(ref_alignment.length());
    while (seq_index < boundary and alignment_index < ref_length) {
        alignment_index++;
        char current_base = ref_alignment[alignment_index];
        if (current_base != '-') { seq_index++;}
    }
    return alignment_index;
}

bool query_contains_ref_bases(Alignment alignment, VarBoundary var_boundary)
{ 
    std::string btop = alignment.btop;
    int64_t ref_start = alignment.ref_start;
    int64_t var_start = var_boundary.start;
    int64_t var_stop = var_boundary.stop;
    // Adjust the variant boundary indices to be compatible with the local alignment given by the BTOP string
    var_start = var_start - ref_start;
    var_stop = var_stop - ref_start;
    // Convert the BTOP string into an alignment
    std::vector<std::string> btop_list = get_btop_list(btop);
    std::string ref_alignment = find_reference_alignment(btop_list);
    // Translate the variant sequence boundaries to an alignment boundaries
    int64_t alignment_start = get_alignment_var_boundary(ref_alignment,var_start);
    int64_t alignment_stop = get_alignment_var_boundary(ref_alignment,var_stop);
    // If any base in the interval containing the variant is not a match, then the query does not contain the variant
    for (int i = alignment_start; i <= alignment_stop; i++) { 
        if (ref_alignment[i] != MATCH_g) {return false;} 
    }
    return true; 
}
