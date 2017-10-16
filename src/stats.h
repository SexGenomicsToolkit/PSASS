#pragma once
#include "variant_data.h"

void analysis(std::ifstream& input_file, std::ofstream& fst_output_file, std::ofstream& snps_output_file, uint16_t min_reads_sex, float range, float min_fst);
