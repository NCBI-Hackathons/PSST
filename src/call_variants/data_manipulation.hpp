#ifndef DATA_MAIPULATION_H
#define DATA_MAIPULATION_H

#include "structs.hpp"

std::vector<std::string> split(const std::string &str);

std::vector< std::vector<BlastOutput> > partition_mbo_paths(const std::vector<BlastOutput> &mbo_paths, int64_t n);

#endif /* DATA_MAIPULATION_H */
