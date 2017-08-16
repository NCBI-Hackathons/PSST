#ifndef DATA_MAIPULATION_H
#define DATA_MAIPULATION_H

std::vector<std::string> split(const std::string &str);

struct Alignment
{
    std::string btop;
    std::string ref_start;
    std::string ref_stop;
};

struct FlankInfo
{
    int64_t start;
    int64_t stop;
};

#endif /* DATA_MAIPULATION_H */
