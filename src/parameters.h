#pragma once
#include <string>
#include <fstream>

struct Parameters{

    uint min_depth;
    float min_fst;
    float snp_range;
    float fixed_range;
    uint window_size;
    uint output_resolution;

    uint male_pool;
    bool output_coverage;
    bool output_snps_pos;
    bool output_genes = false;

    std::string input_file_path;
    std::string gff_file_path;
    std::string output_file_path;

    std::ifstream input_file;
    std::ifstream gff_file;

    std::ofstream snps_output_file;
    std::ofstream snps_pos_output_file;
    std::ofstream fst_threshold_output_file;
    std::ofstream fst_window_output_file;
    std::ofstream coverage_output_file;
    std::ofstream genes_output_file;
};
