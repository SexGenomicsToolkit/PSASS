#pragma once
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <unordered_map>
#include "input_data.h"
#include "gff_file.h"
#include "logs.h"
#include "parameters.h"
#include "pool_data.h"


struct OutputFile {

    std::string suffix = "";
    std::string header = "";
    std::string path = "";
    std::ofstream file;

    OutputFile() {}
    OutputFile(const std::string& suffix, const std::string& header) {this->suffix = suffix; this->header = header;}
    void open(const std::string& prefix, Logs* logs);
};


class OutputHandler {

    public:

        OutputHandler() {}
        OutputHandler(Parameters* parameters, InputData* input_data, PoolBaseData* male_pool, PoolBaseData* female_pool, bool male_index, bool female_index,
                      std::map<std::string, std::map<uint, float[3]>>* depth, std::unordered_map<std::string, Gene>* genes, Logs* logs);
        void output_fst_position(float fst);
        void output_fst_window(float fst_parts[2]);
        void output_snp_position(std::string sex);
        void output_snp_window(uint32_t snps_total[2]);
        void output_depth(float* average_depth);
        void output_genes(float* average_depth);

    private:

        OutputFile fst_position_output_file {"fst_position", "Contig\tPosition\tFst\n"};
        OutputFile fst_window_output_file {"fst_window", "Contig\tPosition\tFst\n"};
        OutputFile snps_position_output_file {"snps_position", "Contig\tPosition\tSex\tM_A\tM_T\tM_G\tM_C\tM_I\tF_A\tF_T\tF_G\tF_C\tF_I\n"};
        OutputFile snps_window_output_file {"snps_window", "Contig\tPosition\tMales\tFemales\n"};
        OutputFile depth_output_file {"depth", "Contig\tPosition\tMales_depth_abs\tFemales_depth_abs\tMales_depth_rel\tFemales_depth_rel\n"};
        OutputFile genes_output_file {"genes", "Contig\tStart\tEnd\tID\tName\tProduct\t"
                                      "M_depth\tM_depth_corr\tM_depth_coding\tM_depth_coding_corr\tM_depth_noncoding\tM_depth_noncoding_corr\t"
                                      "F_depth\tF_depth_corr\tF_depth_coding\tF_depth_coding_corr\tF_depth_noncoding\tF_depth_noncoding_corr\t"
                                      "M_snps\tM_snps_coding\tM_snps_noncoding\tF_snps\tF_snps_coding\tF_snps_noncoding\n"};

        Logs* logs;
        Parameters* parameters;
        InputData* input_data;
        PoolBaseData* male_pool;
        PoolBaseData* female_pool;
        bool male_index = 1;
        bool female_index = 0;
        std::map<std::string, std::map<uint, float[3]>>* depth;
        std::unordered_map<std::string, Gene>* genes;

        void create_output_files();
        void open_output_file(OutputFile& output_file);
};

