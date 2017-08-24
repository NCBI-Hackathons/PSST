#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <map>
#include <cstdlib> // For std::exit
#include <unistd.h> // For getopt
#include <algorithm> // For std::remove
// For multiprocessing
#include <future>
#include <thread>
// Custom libraries
#include "structs.hpp"
#include "data_manipulation.hpp"
#include "queries_with_ref_bases.hpp"
#include "variant_calling.hpp"

void remove_new_lines(std::string &str)
/* Removes new line characters from the given string. */
{
    str.erase( std::remove(str.begin(), str.end(), '\n'), str.end() );
}

std::vector<BlastOutput> get_mbo_paths(std::string mbo_list_path)
/* Retrieves the paths to the MBO alignment files */
{
    // Open the file
    std::ifstream mbo_list_file (mbo_list_path,std::ios::in);
    if ( not mbo_list_file.is_open() ) {
        std::cerr << "Error: unable to read SRA list file.\n";
        std::exit(1);
    }

    std::string line;
    std::vector<BlastOutput> blast_output_list;

    // Each line of the file corresponds to an line number
    while ( std::getline(mbo_list_file,line) ) {
        // Remove new lines from the line
        remove_new_lines(line);
        // Split into tokens
        std::vector<std::string> tokens = split(line); 
        // The first token is the accession, the second is the MBO path
        BlastOutput blast_output;
        blast_output.accession = tokens.at(0);
        blast_output.path = tokens.at(1); 
        // Add it to the list
        blast_output_list.push_back(blast_output); 
    }

    mbo_list_file.close();
    return blast_output_list;
}

std::map<std::string,VarBoundary> get_var_boundary_map(std::string var_boundary_path)
/* Retrieve the variant flanking sequence information from a file */
{
    // Open the file
    std::ifstream var_boundary_file (var_boundary_path,std::ios::in);
    if ( not var_boundary_file.is_open() )
    {
        std::cerr << "Error: unable to read variant information file.\n";
        std::exit(1);
    }

    std::string line;
    std::map<std::string, VarBoundary> var_boundary_map;

    // Each line corresponds to the information for a single variant
    while ( std::getline(var_boundary_file,line) ) {
        // Remove new lines from the line
        remove_new_lines(line);
        // Split into tokens
        std::vector<std::string> tokens = split(line);
        // Put the tokens into the variant info object
        std::string accession = tokens.at(0);
        VarBoundary var_boundary;
        // Convert the strings into integers
        var_boundary.start = std::stoi( tokens.at(1) );
        var_boundary.stop = std::stoi( tokens.at(2) ); 
        var_boundary_map[accession] = var_boundary;
    } 

    var_boundary_file.close();
    return var_boundary_map;
}

SraAlignments get_sra_alignments(BlastOutput blast_output)
/* Retrieves the Magic-BLAST alignments for a particular SRA dataset */
{
    std::string accession = blast_output.accession;
    std::string path = blast_output.path;

    std::ifstream mbo_file (path,std::ios::in);
    if ( not mbo_file.is_open() ) {
        std::cerr << "Error: unable to read " << path << std::endl;
        std::exit(1);
    }

    std::string line;
    std::vector<Alignment> alignments;
    Alignment alignment;
    
    while ( std::getline(mbo_file,line) ) {
        // Remove new lines from the line
        remove_new_lines(line);
        // Split into tokens like awk
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
        alignment.variant_acc = variant_acc;
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

void write_tsv(const std::vector<CalledVariants> &called_variants_list, std::string output_path)
/* Create a TSV file describing which variants exist in each SRA dataset */
{
    std::ofstream tsv_file (output_path, std::ios::out | std::ios::trunc);
    if (tsv_file.is_open()) {
        tsv_file << "SRA Dataset \t Heterozygous \t Homozygous\n";
        for (size_t i = 0; i < called_variants_list.size(); i++) {
            CalledVariants sra_variants = called_variants_list.at(i);
            std::string accession = sra_variants.accession;
            std::vector<std::string> heterozygous = sra_variants.heterozygous;
            std::vector<std::string> homozygous = sra_variants.homozygous;
            tsv_file << accession << " \t ";
            for (size_t het_index = 0; het_index < heterozygous.size(); het_index++) {
                tsv_file << heterozygous.at(het_index) << " ";
            }
            tsv_file << " \t ";
            for (size_t hom_index = 0; hom_index < homozygous.size(); hom_index++) {
                tsv_file << homozygous.at(hom_index) << " ";
            } 
            tsv_file << "\n";
        }
        tsv_file.close();
    } else {
        std::cerr << "Error: unable to create TSV file.\n";
    }
} 

int main(int argc, char *argv[])
{
    int opt;
    std::string help = "Description: Given a list of paths to Magic-BLAST output files and a file containing information concerning variant boundaries, returns a TSV file describing which heterozygous and homozygous variants exist in each dataset.\n"; 
    std::string usage = "[-h help and usage] [-m MBO list path] [-v var boundary info path] [-o output path]\n";
    

    std::string mbo_list_path = "";
    std::string var_boundary_path = "";
    std::string output_path = "";

    while ( (opt = getopt(argc, argv, "hm:v:o:")) != 1 ) {
        switch (opt) {
            case 'h':
                std::cout << help;
                std::cout << usage;
                return 0;
            case 'm':
                mbo_list_path = optarg;
                break;
            case 'v':
                var_boundary_path = optarg;     
                break;
            case 'o':
                output_path = optarg;
                break;
            default:
                std::cerr << "Error: unrecognized option.\n";
                std::cerr << usage;
                return 1;
        }
    }

    bool opts_incomplete = false;

    if (mbo_list_path == "") {
        std::cerr << "Error: please provide the path to the file containing the paths of the MBO files.\n";
        opts_incomplete = true;
    }
    if (var_boundary_path == "") {
        std::cerr << "Error: please provide the path to the variant flanking sequence information file.\n";
        opts_incomplete = true;
    }
    if (output_path == "") {
        std::cerr << "Error: please specify an output path.\n";
        opts_incomplete = true;
    }
    if (opts_incomplete) {
        std::cerr << usage;
        return 1;
    }

    std::vector<BlastOutput> mbo_paths_list = get_mbo_paths(mbo_list_path);    
    std::map<std::string,VarBoundary> var_boundary_map = get_var_boundary_map(var_boundary_path);
    std::vector<SraAlignments> sra_alignments_list;
    // Get the alignments for each SRA dataset
    for (size_t i = 0; i < mbo_paths_list.size(); i++) {
        BlastOutput blast_output = mbo_paths_list.at(i);
        SraAlignments sra_alignments = get_sra_alignments(blast_output); 
        sra_alignments_list.push_back(sra_alignments);
    }
    // Call them variants
    std::vector<CalledVariants> called_variants_list = call_variants(sra_alignments_list,var_boundary_map); 
    // Write the output file
    write_tsv(called_variants_list,output_path);
}
