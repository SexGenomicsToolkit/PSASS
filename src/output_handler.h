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


struct OutputFile {

    std::string suffix = "";
    std::string header = "";
    std::string path = "";
    std::ofstream file;

    OutputFile() {}
    OutputFile(const std::string& suffix) {this->suffix = suffix;}
    void open(const std::string& prefix);
};


class OutputHandler {

    public:

        OutputHandler() {}
        OutputHandler(Parameters* parameters, InputData* input_data, PairBaseData* pair_data, std::map<std::string, std::map<uint, float[3]>>* depth, std::unordered_map<std::string, Gene>* genes);
        void output_fst_position(float fst);
        void output_fst_window(float fst_parts[2]);
        void output_snp_position(std::string& pool_id);
        void output_snp_window(uint32_t snps_total[2]);
        void output_depth(float* average_depth);
        void output_genes(float* average_depth);

    private:

        OutputFile fst_position_output_file {"fst_position"};
        OutputFile fst_window_output_file {"fst_window"};
        OutputFile snps_position_output_file {"snps_position"};
        OutputFile snps_window_output_file {"snps_window"};
        OutputFile depth_output_file {"depth"};
        OutputFile genes_output_file {"genes"};

        Parameters* parameters;
        InputData* input_data;
        PairBaseData* pair_data;
        std::map<std::string, std::map<uint, float[3]>>* depth;
        std::unordered_map<std::string, Gene>* genes;

        void create_output_files();
        void open_output_file(OutputFile& output_file);
};

