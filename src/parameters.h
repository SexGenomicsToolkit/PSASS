#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "utils.h"

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

    // Pileup options
    std::vector<std::string> alignment_files;
    std::string reference_file = "";
    uint min_mapping_quality = 0;

    // Kpool options
    std::string table1_file_path = "";
    std::string table2_file_path = "";
    std::string tmp_file_prefix = "";
    std::string output_prefix = "";
    uint32_t min_kmer_presence_depth = 0;
    uint32_t max_kmer_presence_depth = 99999;
    uint32_t max_kmer_absence_depth = 0;
    std::string pool_to_filter = "";

    void print_psass() {

        log("Parameter values:");
        std::cerr << "Input:\n";
        std::cerr << "    - Input file:     " << "" << this->input_file_path << "\n";
        std::cerr << "    - Input format:   ";
        this->popoolation_format ? std::cerr <<  "popoolation2\n" : std::cerr << "psass\n";
        std::cerr << "    - GFF file:       ";
        this->gff_file_path != "" ?  std::cerr << "\"" << this->gff_file_path << "\"\n" : std::cerr << "None\n";
        std::cerr << "    - Pool 1 name:    " << this->pool1_id << "\n";
        std::cerr << "    - Pool 2 name:    " << this->pool2_id << "\n";
        std::cerr << "Output:\n";
        std::cerr << "    - Size of sliding window:               " << this->window_size << "\n";
        std::cerr << "    - Output resolution:                    " << this->output_resolution << "\n";
        std::cerr << "    - Output file for sliding window:       " << this->output_file_path << "\n";
        std::cerr << "    - Output file for pool-specific SNPs:   ";
        this->snp_pos_file_path != "" ? std::cerr << this->snp_pos_file_path << "\n" : std::cerr << "None\n";
        std::cerr << "    - Output file for high FST positions:   ";
        this->fst_pos_file_path != "" ? std::cerr << this->fst_pos_file_path << "\n" : std::cerr << "None\n";
        std::cerr << "    - Output file for gene metrics:   ";
        this->genes_file_path != "" ? std::cerr << this->genes_file_path << "\n" : std::cerr << "None\n";
        std::cerr << "Analysis:\n";
        std::cerr << "    - Minimum depth to consider a site:                    " << this->min_depth << "\n";
        std::cerr << "    - Minimum fst to output high fst site:                 " << this->min_fst << "\n";
        std::cerr << "    - Allele frequency of a heterozygous locus:            " << this->freq_het << "\n";
        std::cerr << "    - Allele frequency of a homozygous locus:              " << this->freq_hom << "\n";
        std::cerr << "    - Range of allele frequency of a heterozygous locus:   " << this->range_het << "\n";
        std::cerr << "    - Range of allele frequency of a homozygous locus:     " << this->range_hom << "\n";
        std::cerr << "    - Count contiguous snps as one: " << std::boolalpha << this->group_snps << "\n";
    }
};
