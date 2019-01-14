#pragma once
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include "input_data.h"
#include "pool_data.h"
#include "parameters.h"


struct OutputFile {

    std::string suffix;
    std::string header;
    std::string path;
    std::ofstream file;

    OutputFile() {}
    OutputFile(const std::string& suffix, const std::string& header) {this->suffix = suffix; this->header = header;}
    void open(const std::string& prefix);
};


class OutputHandler {

    public:

        OutputHandler() {}
        OutputHandler(Parameters* parameters, InputData* input_data, PoolBaseData* male_pool, PoolBaseData* female_pool, bool male_index, bool female_index);
        void output_snp_position(std::string sex);
        void output_snp_window(uint16_t snps_total[2]);

    private:

        OutputFile fst_position_output_file {"fst_position", "Contig\tPosition\tFst\n"};
        OutputFile fst_window_output_file {"fst_window", "Contig\tPosition\tFst\n"};
        OutputFile snps_position_output_file {"snps_position", "Contig\tPosition\tSex\tM_A\tM_T\tM_G\tM_C\tM_I\tF_A\tF_T\tF_G\tF_C\tF_I\n"};
        OutputFile snps_window_output_file {"snps_window", "Contig\tPosition\tMales\tFemales\n"};
        OutputFile coverage_output_file {"coverage", "Contig\tPosition\tMales_rel\tFemales_rel\tMales_abs\tFemales_abs\n"};
        OutputFile genes_output_file {"genes", "Contig\tStart\tEnd\tName\tProduct\t"
                                      "Cov_males\tCov_males_corr\tCov_males_coding\tCov_males_coding_corr\tCov_males_noncoding\tCov_males_noncoding_corr\t"
                                      "Cov_females\tCov_females_corr\tCov_females_coding\tCov_females_coding_corr\tCov_females_noncoding\tCov_females_noncoding_corr\t"
                                      "Snp_males\tSnp_males_coding\tSnp_males_noncoding\tSnp_females\tSnp_females_coding\tSnp_females_noncoding\n"};

        InputData* input_data;
        PoolBaseData* male_pool;
        PoolBaseData* female_pool;
        bool male_index;
        bool female_index;
        Parameters* parameters;

        void create_output_files();
        void open_output_file(OutputFile& output_file);
};

