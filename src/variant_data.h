#pragma once
#include "utils.h"

struct Position {
    uint16_t male_bases[5] {0, 0, 0, 0, 0};
    uint16_t female_bases[5] {0, 0, 0, 0, 0};
};

typedef std::unordered_map<std::string, std::vector<Position>> variants;
