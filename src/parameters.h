#pragma once
#include <string>
#include <fstream>

struct Parameters{

    // Analysis parameters
    uint min_depth = 10;
    float min_fst = float(0.1);
    float freq_het = float(0.5);
    float freq_hom = float(1.0);
    float range_het = float(0.15);
    float range_hom = float(0.05);
    float min_het = float(0.35);
    float max_het = float(0.65);
    float min_hom = float(0.95);
    uint window_size = 100000;
    uint window_range = 50000;
    uint output_resolution = 1000;
    bool group_snps= false;

    // Output options
    bool output_fst_pos = true;
    bool output_fst_win = true;
    bool output_snps_pos = true;
    bool output_snps_win = true;
    bool output_depth = true;
    bool output_genes = false;

    // Input options
    std::string input_file_path = "";
    std::string gff_file_path = "";
    std::string output_prefix = "";
    uint male_pool = 1;

    // Input file streams
    std::ifstream input_file;
    std::ifstream gff_file;
};
