#pragma once
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <unordered_map>
#include "input_data.h"
#include "gff_file.h"
#include "pair_data.h"
#include "parameters.h"
#include "pool_data.h"


class OutputHandler {

    public:

        OutputHandler() {}
        OutputHandler(Parameters& parameters);
        void output_window(std::map<std::string, std::map<uint, float[6]>>& output_data, float* average_depth, std::unordered_map<std::string, uint64_t>& contig_lengths);
        void output_fst(float fst, InputData& input_data);
        void output_snp(std::string& pool_id, PairBaseData& pair_data, InputData& input_data);
        void output_genes(std::unordered_map<std::string, Gene>& genes, float* average_depth);

    private:

        std::ofstream window_output_file;
        std::ofstream fst_position_output_file;
        std::ofstream snp_position_output_file;
        std::ofstream genes_output_file;

        void create_output_files();
        void open_output_file(std::ofstream& output_file, std::string& path);
};

