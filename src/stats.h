#pragma once
#include "variant_data.h"

void calculate_fst(variants &data, std::ofstream& output_file, uint32_t window, std::unordered_map<std::__cxx11::string, uint32_t>& contig_lengths);
