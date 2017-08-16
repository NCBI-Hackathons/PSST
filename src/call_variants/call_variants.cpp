#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cstdlib> // For std::exit
#include <map>
// For multithreading
#include <future>
#include <thread>
// Custom libraries
#include "data_manipulation.hpp"
#include "queries_with_ref_bases.hpp"

struct BlastOutput
{
    std::string accession;
    std::string mbo_path;
};

struct VarInfo
{
    int64_t start;
    int64_t stop;
};

struct VariantFreq
{
    std::string accession;
    int64_t matches;
    int64_t mismatches;
};

struct CalledVariants
{
    std::vector<std::string> homozygous;
    std::vector<std::string> heterozygous;
};

struct SraVariants
{
    std::string accession;
    CalledVariants called_variants;
};

struct SraAlignments
{
    std::string accession;
    std::vector<Alignment> alignments;
};

std::vector< std::string > get_sra_paths(std::string sra_list_path)
// Retrieves the paths to the MBO alignment files
{
    std::ifstream sra_list_file (sra_list_path,std::ios::in);
    if ( not sra_list_file.is_open() )
    {
        std::cerr << "Error: unable to read SRA list file.\n";
        std::exit(1);
    }

    std::string line;
    std::vector<BlastOutput> sra_paths;

    // Each line of the file corresponds to an line number
    while ( std::getline(sra_list_file,line) ) {
        // Remove new lines from the line
        line.erase( std::remove(line.begin(), line.end(), '\n'), line.end() );
        // Split into tokens
        std::vector<std::string> tokens = split(line); 
        // The first token is the accession, the second is the MBO path
        BlastOutput sra_path;
        sra_path.accession = tokens.at(0);
        sra_path.mbo_path = tokens.at(1); 
        // Add it to the list
        sra_paths.push_back(sra_path); 
    }
    return sra_paths
}

std::map<std::string, VarInfo> get_var_info_map(std::string var_info_path)
{
    std::ifstream var_info_file (var_info_path,std::ios::in);
    if ( not var_info_file.is_open() )
    {
        std::cerr << "Error: unable to read variant information file.\n";
        std::exit(1);
    }

    std::string line;
    std::map<std::string, VarInfo> var_info_map;

    // Each line corresponds to the information for a single variant
    while ( std::getline(var_info_file,line) ) {
        // Remove new lines from the line
        line.erase( std::remove(line.begin(), line.end(), '\n'), line.end() );
        // Split into tokens
        std::vector<std::string> tokens = split(line);
        // Put the tokens into the variant info object
        accession = tokens.at(0);
        VarInfo var_info;
        // Convert the strings into integers
        var_info.start = std::stoi( tokens.at(1) );
        var_info.stop = std::stoi( tokens.at(2) ); 
        var_info_map[accession] = var_info;
    } 

    return var_info_map;
}

SraAlignments get_sra_alignments(BlastOutput blast_output)
{
    accession = blast_output.accession
    mbo_path = blast_output.mbo_path;

    std::ifstream mbo_file (mbo_path,std::ios::in);
    if ( not mbo_file.is_open() ) {
        std::cerr << "Error: unable to read " << mbo_path << std::endl;
        std::exit(1);
    }

    std::string line;
    std::vector<Alignment> alignments;
    Alignment alignment;
    
    while ( std::getline(mbo_file,line) ) {
        // Remove new lines from the line
        line.erase( std::remove(line.begin(), line.end(), '\n'), line.end() );
        // Split into tokens
        std::vector<std::string> tokens = split(line);
        std::string variant_acc = tokens.at(1);
        int ref_start = std::stoi( tokens.at(8) );
        int ref_stop = std::stoi( tokens.at(9) );
        if (ref_start > ref_stop) {
            int temp = ref_start;
            ref_start = ref_stop;
            ref_stop = temp;
        }
        std::string btop = tokens.at(16);
        alignment.ref_start = ref_start;
        alignment.ref_stop = ref_stop;
        alignment.btop = btop;
        alignments.push_back(alignment); 
    } 
    SraAlignments sra_alignments;
    sra_alignments.accession = accession;
    sra_alignments.alignments = alignments;
    return sra_alignments;
} 

CalledVariants call_variants(std::vector<VariantFreq> var_freq_list)
{
    CalledVariants called_variants;
    std::vector<std::string> homozygous;
    std::vector<std::string> heterozygous; 
    for (int64_t index = 0; index < var_freq.size(); index ++) {
        var_freq = var_freq_list.at(index);
        accession = var_freq.accession;
        matches = var_freq.matches;
        mismatches = var_freq.mismatches;
        total_reads = matches + mismatches;
        if (total_reads != 0) {
            float percent_match = matches/total_reads;
            if (percent_match > 0.8) {
                homozygous.push_back(accession);
            } else if (percent_match > 0.3) {
                heterozygous.push_back(accession);
            }
        }
    }
    called_variants.homozygous = homozygous;
    called_variants.heterozygous = heterozygous;
    return called_variants;
}

std::vector<SraVariants> call_variants_in_datasets(std::vector<BlastOutput> blast_outputs, 
                                                   std::map<std::string,VarInfo> var_info_map)
{
    for (int index = 0; index < blast_outputs.size(); index++) {
        BlastOutput blast_output = blast_outputs.get(index);
        SraAlignments sra_alignments = get_sra_alignments(blast_output);
        std::map<std::string,
    } 
} 
