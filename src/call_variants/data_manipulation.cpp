#include <string>
#include <sstream>
#include <vector>
#include "structs.hpp"
#include "data_manipulation.hpp"

std::vector<std::string> split(const std::string &str)
{
    std::vector<std::string> elems;
    std::stringstream ss(str);
    std::string item;

    while (std::getline(ss, item, ' ')) {
        if (!item.empty()) {elems.push_back(item);}
    }
    return elems;
}

std::vector< std::vector<BlastOutput> > partition_mbo_paths(std::vector<BlastOutput> &mbo_paths, int64_t n)
/* Partitions the mbo_paths vector into n partitions of near identical size */
{
    std::vector<BlastOutput> partition;
    std::vector< std::vector<BlastOutput> > partitions;
    std::vector<BlastOutput>::iterator iter = mbo_paths.begin();
    int64_t partition_size = mbo_paths.size() / n;
    // This is the number of partitions that will have size partition_size + 1
    int64_t extra = mbo_paths.size() % n;
    // Partition the vector
    for (int i = 0; i < n; i++) {
        int64_t offset = partition_size + (i < extra ? 1 : 0);
        partition.assign(iter,iter+offset);
        partitions.push_back(partition);
        iter += offset;
    }
    return partitions;
}
