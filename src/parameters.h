#pragma once
#include <string>
#include <fstream>

struct Parameters{

    // Analysis parameters
    uint min_depth;
    float min_fst;
    float freq_het;
    float freq_hom;
    float range_het;
    float range_hom;
    float min_het;
    float max_het;
    float min_hom;
    uint window_size;
    uint window_range;
    uint output_resolution;
    bool ignore_indels;

    // Output options
    bool output_fst_pos;
    bool output_fst_win;
    bool output_snps_pos;
    bool output_snps_win;
    bool output_depth;
    bool output_genes = false;

    // Input options
    std::string input_file_path;
    std::string gff_file_path;
    std::string output_prefix;
    uint male_pool;

    // Input file streams
    std::ifstream input_file;
    std::ifstream gff_file;

};
