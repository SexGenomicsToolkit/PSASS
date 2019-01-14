#pragma once
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>
#include "input_data.h"
#include "pool_data.h"
#include "parameters.h"


struct OutputFile {

    std::string suffix = "";
    std::string header = "";
    std::string path = "";
    std::ofstream file;

    OutputFile() {}
    OutputFile(const std::string& suffix, const std::string& header) {this->suffix = suffix; this->header = header;}
    void open(const std::string& prefix);
};


class OutputHandler {

    public:

        OutputHandler() {}
        OutputHandler(Parameters* parameters, InputData* input_data, PoolBaseData* male_pool, PoolBaseData* female_pool, bool male_index, bool female_index);
        void output_fst_position(float fst);
        void output_snp_position(std::string sex);
        void output_snp_window(uint16_t snps_total[2]);
        void output_depth(std::map<std::string, std::map<uint, float[2]>>& depth, uint64_t* total_depth, uint64_t& total_bases);

    private:

        OutputFile fst_position_output_file {"fst_position", "Contig\tPosition\tFst\n"};
        OutputFile fst_window_output_file {"fst_window", "Contig\tPosition\tFst\n"};
        OutputFile snps_position_output_file {"snps_position", "Contig\tPosition\tSex\tM_A\tM_T\tM_G\tM_C\tM_I\tF_A\tF_T\tF_G\tF_C\tF_I\n"};
        OutputFile snps_window_output_file {"snps_window", "Contig\tPosition\tMales\tFemales\n"};
        OutputFile depth_output_file {"depth", "Contig\tPosition\tMales_depth_abs\tFemales_depth_abs\tMales_depth_rel\tFemales_depth_rel\n"};
        OutputFile genes_output_file {"genes", "Contig\tStart\tEnd\tName\tProduct\t"
                                      "Males_depth\tMales_depth_corr\tMales_depth_coding\tMales_depth_coding_corr\tMales_depth_noncoding\tMales_depth_noncoding_corr\t"
                                      "Females_depth\tFemales_depth_corr\tFemales_depth_coding\tFemales_depth_coding_corr\tFemales_depth_noncoding\tFemales_depth_noncoding_corr\t"
                                      "Males_snps\tMales_snps_depth_coding\tMales_snps_depth_noncoding\tFemales_snps\tFemales_snps_depth_coding\tFemales_snps_depth_noncoding\n"};

        InputData* input_data;
        PoolBaseData* male_pool;
        PoolBaseData* female_pool;
        bool male_index = 1;
        bool female_index = 0;
        Parameters* parameters;

        void create_output_files();
        void open_output_file(OutputFile& output_file);
};

