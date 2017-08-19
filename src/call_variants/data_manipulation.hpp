#ifndef DATA_MAIPULATION_H
#define DATA_MAIPULATION_H

std::vector<std::string> split(const std::string &str);

struct Alignment
{
    std::string variant_acc;
    std::string btop;
    int64_t ref_start;
    int64_t ref_stop;
};

struct VarInfo
{
    int64_t start;
    int64_t stop;
};

#endif /* DATA_MAIPULATION_H */
