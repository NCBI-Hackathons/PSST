#include <string>
#include <algorithm>
#include <cctype>
#include "queries_with_ref_bases.hpp"
#include "data_manipulation.hpp"

// Global variables
std::string GAP_g = "_";
std::string MATCH_g = "=";
int REF_g = 1;

std::vector<std::string> get_btop_list(std::string btop) 
{   // We treat introns as gaps
    std::replace( btop.begin(), btop.end(), '^', '_' );
    std::string temp_btop = "";
    int num_base = 0;
    bool in_gap = false;
    bool prev_was_int = false;
    for (int i = 0; i < btop.length(); i++) {
        std::string current_char = btop[i];
        if (std::isalpha(current_char) or current_char == "-") {
            num_base ++;
            if (num_base % 2 == 1) { temp_btop = temp_btop + " "; }
            prev_was_int = false;
        } else if (current_char == "_") {
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
/* Given a BTOP list, returns a pseudoalignment of the reference */
{ 
    std::string ref_alignment = "";
    std::string base; 
    bool gap = false;
    for (int i = 0; i < btop_list.size(); i++) {
        std::string op = btop_list.at(i);
        if (op.find("_") != std::string::npos) {
            gap = true;
            op.erase(std::remove(op.begin(),op.end(),"_"),op.end());
        } else {
            gap = false;
        }
        if (std::isdigit(op)) {
            int num_matches = std::stoi(btop);
            if (gap) { base = GAP_g; } else { base = MATCH_g; }
            for (int u = 0; u < num_matches; u++) { ref_alignment = ref_alignment + base; }
        } else { ref_alignment = ref_alignment + op[REF_g]; }
    }
    return ref_alignment;
}

int64_t translate_var_boundary(std::string ref_alignment, int64_t boundary) 
/* Translates modifies the variant boundary to be compatible with the reference pseudoalignment */
{ 
    int alignment_boundary = 0;
    int seq_index = 0;
    while (seq_index < boundary and alignment_boundary < ref_alignment.length()) {
        std::string current_base = ref_alignment[alignment_boundary];
        if (current_base != "-") { seq_index++; }
        alignment_boundary++;
    }
    return alignment_boundary;
}

bool query_contains_ref_bases(AlignmentInfo alignment_info, FlankInfo flank_info)
/* Determines whether the given read contains the given variant */
{ 
    std::string btop = alignment_info.btop;
    int64_t ref_start = alignment.ref_start;
    int64_t ref_stop = alignment.ref_stop;
    int64_t var_start = alignment.ref_start;
    int64_t var_stop = alignment.ref_stop;
    // Adjust the variant boundary indices to be compatible with the local alignment given by the BTOP string
    var_start = var_start - ref_start;
    var_stop = var_stop - ref_stop;
    // Convert the BTOP string into an alignment
    std::vector<std::string> btop_list = get_btop_list(btop);
    std::string ref_alignment = find_reference_alignment(btop_list);
    // Translate the variant sequence boundaries to an alignment boundaries
    int64_t alignment_start = translate_var_boundary(var_start);
    int64_t alignment_stop = translate_var_boundary(var_stop);
    // If any base in the interval containing the variant is not a match, then the query does not contain the variant
    for (int i = alignment_start; i < alignment_stop; i++) { if (ref_alignment[i] != MATCH_g) {return false;} }
    return true; 
}
