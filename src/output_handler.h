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


struct PointOutputData {

    PairBaseData base_data;
    std::string contig;
    uint64_t position;
    uint8_t type = 0;  // 0 -> fst, 1 --> snp in pool1, 2 --> snp in pool2

    PointOutputData(PairBaseData base_data, std::string contig, uint64_t position, uint8_t type) {
        this->base_data = base_data;
        this->contig = contig;
        this->position = position;
        this->type = type;
    }
};


class OutputHandler {

    public:

        OutputHandler() {}
        OutputHandler(Parameters& parameters);
        void output_window(std::map<std::string, std::map<uint, float[6]>>& output_data, float* average_depth, std::unordered_map<std::string, uint64_t>& contig_lengths);
        void output_bases(std::vector<PointOutputData> base_output_data, std::unordered_map<std::string, uint64_t>& contig_lengths);
        void output_genes(std::unordered_map<std::string, Gene>& genes, float* average_depth);

    private:

        std::ofstream window_output_file;
        std::ofstream fst_position_output_file;
        std::ofstream snp_position_output_file;
        std::ofstream genes_output_file;

        std::string pool_id[2];
        uint min_depth = 0;

        void create_output_files();
        void open_output_file(std::ofstream& output_file, std::string& path);
};

