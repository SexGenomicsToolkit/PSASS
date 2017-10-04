#pragma once
#include "variant_data.h"

void get_contig_lengths(std::ifstream& vcf_file, std::unordered_map<std::string, uint32_t>& contig_lengths);
void get_variant_data(std::ifstream& file, variants& data, uint16_t window, std::unordered_map<std::__cxx11::string, uint32_t> &contig_lengths);

