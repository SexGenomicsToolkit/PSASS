#pragma once
#include <string>
#include <fstream>

struct Parameters{

    std::string command = "";

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
    uint output_resolution = 10000;
    bool group_snps= false;

    // Output options
    std::string output_file_path = "";
    std::string snp_pos_file_path = "";
    std::string fst_pos_file_path = "";
    std::string genes_file_path = "";

    // Input options
    std::string input_file_path = "";
    std::string gff_file_path = "";
    std::string pool1_id = "females";
    std::string pool2_id = "males";
    bool popoolation_format = false;
};
