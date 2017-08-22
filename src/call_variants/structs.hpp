#ifndef STRUCTS_H
#define STRUCTS_H

#include <map>

struct Alignment
/* Holds alignment information outputted by Magic-BLAST */
{
    std::string variant_acc;
    std::string btop;
    int64_t ref_start;
    int64_t ref_stop;
};

struct VarBoundary
/* Specifies where the variant lies in the variant flanking sequence */
{
    int64_t start;
    int64_t stop;
};

struct BlastOutput
/* Information concerning where to find the magicblast output path for a particular SRA dataset */
{
    std::string accession;
    std::string path;
};

struct CalledVariants
/* Holds lists of the homozygous and heterozygous variants called in a particular SRA dataset */
{
    std::string accession;
    std::vector<std::string> homozygous;
    std::vector<std::string> heterozygous;
};

struct SraAlignments
/* Container for the Magic-BLAST alignments for a particular SRA dataset */
{
    std::string accession;
    std::vector<Alignment> alignments;
};

struct VariantFrequency
/* Container for the number of matched and mismatched reads that span a particular variant */
{
    int64_t matches;
    int64_t mismatches;
};

struct VariantFrequenciesMap
/* The set of all variants that an SRA dataset is mapped to and their identity frequencies */
{
    std::vector<std::string> variant_accessions;
    std::map<std::string,VariantFrequency> map;
};

#endif /* STRUCTS_H */
