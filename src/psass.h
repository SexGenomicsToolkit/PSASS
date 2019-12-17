#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <iterator>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include "arg_parser.h"
#include "gff_file.h"
#include "input_data.h"
#include "pair_data.h"
#include "parameters.h"
#include "output_handler.h"


// Single base information for a pair of pools (sliding window calculations)
struct WindowBaseData {

    bool snps[2] {0, 0};  // Pair of boolean for sex-specific snps in (pool1, pool2)
    uint32_t depth[2] {0, 0};  // Pair of int for coverage in (pool1, pool2)
    float fst_parts[2] {0, 0};
};


// Sliding window data for a pair of pools
struct Window {

    std::deque<WindowBaseData> data;  // Sliding window data as a deque
    uint32_t snps_in_window[2] {0, 0};  // Pair of int for average snps in current window for (pool1, pool2)
    uint64_t depth_in_window[2] {0, 0};  // Pair of int for average coverage in current window for (pool1, pool2)
    float fst_parts[2] {0, 0};
};


class Psass {

    public:

        std::chrono::steady_clock::time_point t_begin; // Starting time to compute total runtime

        Parameters parameters;  // Parameters updated by the arguments parser

        GFFData gff_data;

        InputData input_data;  // Data related to input file parsing
        OutputHandler output_handler;  // Object handling all output functions

        PairBaseData pair_data;  // PairBaseData object containing information about each pool as well as Fst for this base
        WindowBaseData window_base_data;  // Object containing information about a single base in the sliding window
        Window window;  // Sliding window object

        std::map<std::string, std::map<uint, float[3]>> depth_data;  // Coverage per base for entire genome (needed for relative coverage)
        uint64_t total_depth[2] {0, 0};  // Total coverage (needed for relative coverage)
        uint64_t total_bases = 0;  // Total bases (needed for relative coverage)
        float average_depth[2] = {0.0, 0.0};  // Average depth in male and female pool

        bool consecutive_snps[2] {false, false};

        Psass(int argc, char *argv[]);
        void count_lines();
        void update_nucleotides();
        void update_fst_parts();
        void update_depth();
        void update_snps();
        void update_window(bool end = false);
        void update_genes();
        void output_window_step();
        void process_contig_end();
        void process_line();
        void process_popoolation_field();
        void process_popoolation_subfield();
        void process_psass_field();
        void run();
};








