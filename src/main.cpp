#include <iostream>
#include "variant_data.h"
#include "vcf_parsing.h"
#include "output.h"
#include "stats.h"


int main(int argc, char *argv[]) {

    if (argc != 6) {
        std::cout << "Usage: poolsex sync_file output_file min_reads range min_fst" << std::endl;
        exit(0);
    }

    std::string input_file_path = argv[1];
    std::ifstream input_file;
    input_file.open(input_file_path);

    std::string output_file_path = argv[2];
    std::string fst_output_file_path = output_file_path + "_fst.tsv";
    std::string snps_output_file_path = output_file_path + "_snps.tsv";
    std::ofstream fst_output_file, snps_output_file;
    fst_output_file.open(fst_output_file_path);
    snps_output_file.open(snps_output_file_path);

    uint16_t min_reads_sex = std::stoi(argv[3]);
    float range = std::stof(argv[4]);
    float min_fst = std::stof(argv[5]);

    if (not input_file.is_open()) {

        std::cout << "Error: cannot open input file." << std::endl;
        exit(0);
    }

    if (not fst_output_file.is_open()) {

        std::cout << "Error: cannot open Fst output file." << std::endl;
        exit(0);
    }

    if (not snps_output_file.is_open()) {

        std::cout << "Error: cannot open Snps output file." << std::endl;
        exit(0);
    }

    std::cout << "Analyzing" << std::endl;
    analysis(input_file, fst_output_file, snps_output_file, min_reads_sex, range, min_fst);

    return 0;
}
