#pragma once
#include <stdint.h>
#include <numeric>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <deque>
#include <map>


class PoolBaseData {

    public:

        uint16_t nucleotides[6] {0, 0, 0, 0, 0, 0};  // Raw nucleotides count for a position
        float frequencies[6] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  // Nucleotides frequencies for a position
        float pi = 0.0;  // Pi-statistic
        uint16_t depth = 0;  // Total coverage at this position

        PoolBaseData();

        void compute_total();
        void compute_frequencies();
        void compute_pi();
        void update();

        std::string output_nucleotides();
        std::string output_frequencies();
};
