#pragma once
#include "utils.h"

struct Position {
    uint16_t male_coverage = 0;
    uint16_t female_coverage = 0;
    uint16_t male_heterozygote = 0;
    uint16_t female_heterozygote = 0;
};

typedef std::unordered_map<std::string, std::vector<Position>> variants;
