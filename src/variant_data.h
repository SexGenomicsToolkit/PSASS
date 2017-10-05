#pragma once
#include "utils.h"

struct Position {
    uint16_t male_coverage = 0;
    uint16_t female_coverage = 0;
    bool male_heterozygote = 0;
    bool female_heterozygote = 0;
    bool male_specific_heterozygote = 0;
    bool female_specific_heterozygote = 0;
    uint16_t male_bases[4] {0, 0, 0, 0};
    uint16_t female_bases[4] {0, 0, 0, 0};
    double fst = 0;
    double fisher_p = 0;
    double cochran_p = 0;
};

typedef std::unordered_map<std::string, std::vector<Position>> variants;
