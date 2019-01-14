#pragma once
#include "utils.h"

struct Gene {

    std::string contig;
    std::string start;
    std::string end;
    std::string name;
    std::string product;
    uint coding_length;
    uint noncoding_length;
    uint coverage[6]; // Coding male, Non-coding male, Coding female, Non-coding female, Total male, Total female
    uint snps[6]; // Coding male, Non-coding male, Coding female, Non-coding female, Total male, Total female
};

void read_gff_file(std::ifstream& input_file, std::unordered_map<std::string, std::vector<std::vector<std::string>>>& file, std::unordered_map<std::string, Gene>& genes);
