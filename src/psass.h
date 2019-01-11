#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include <deque>
#include <algorithm>
#include <numeric>
#include <map>
#include <iomanip>
#include "parameters.h"
#include "gff_file.h"
#include "utils.h"
#include "pair_data.h"
#include "arg_parser.h"
#include "analysis.h"
#include "psass.h"


// Single base information for a pair of pools (sliding window calculations)
struct WindowBaseData {

    bool snps[2];  // Pair of boolean for sex-specific snps in (pool1, pool2)
    uint16_t coverage[2];  // Pair of int for coverage in (pool1, pool2)
};


// Sliding window data for a pair of pools
struct Window {

    std::deque<WindowBaseData> data;  // Sliding window data as a deque
    uint16_t snps_average[2];  // Pair of int for average snps in current window for (pool1, pool2)
    uint16_t coverage_average[2];  // Pair of int for average coverage in current window for (pool1, pool2)
};


class Psass {

    public:

        Parameters parameters;

        PoolBaseData* male_pool;
        PoolBaseData* female_pool;
        bool male_index;
        bool female_index;

        std::unordered_map<std::string, std::vector<std::vector<std::string>>> gff_data;
        std::unordered_map<std::string, Gene> genes;
        std::unordered_map<uint, std::pair<std::string, bool>> regions;

        PairBaseData pair_data;
        WindowBaseData window_base_data;
        Window window;

        std::map<std::string, std::map<uint, float[2]>> coverage;  // Coverage per base for entire genome (needed for relative coverage)
        uint64_t total_coverage;  // Total coverage (needed for relative coverage)
        uint64_t total_bases;  // Total bases (needed for relative coverage)

        Psass(int argc, char *argv[]);
        void update_snps();
        void update_window();
        void update_coverage();
        void output_snp(std::string sex, std::string& contig, uint position);
        void run();
};








