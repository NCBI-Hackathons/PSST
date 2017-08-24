#ifndef DATA_MAIPULATION_H
#define DATA_MAIPULATION_H

#include "structs.hpp"

std::vector<std::string> split(const std::string &str);
/* Returns a vector of the words in str where ' ' is the separator */

std::vector< std::vector<BlastOutput> > partition_mbo_paths(const std::vector<BlastOutput> &mbo_paths, int64_t n);
/* Partitions the the vector of BlastOutput objects mbo_paths into n nearly identically-sized vectors, all of which
 * are stored in another vector. */

#endif /* DATA_MAIPULATION_H */
