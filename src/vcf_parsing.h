#pragma once
#include "variant_data.h"

void get_contig_lengths(std::ifstream& vcf_file, std::unordered_map<std::string, uint32_t>& contig_lengths);
void get_variant_data(std::ifstream& male_file, std::ifstream& female_file, std::output& output_file, uint16_t window);

