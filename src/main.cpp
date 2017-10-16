#include <iostream>
#include "variant_data.h"
#include "vcf_parsing.h"
#include "output.h"
#include "stats.h"


int main(int argc, char *argv[]) {

    if (argc < 6) {
        std::cout << "Usage: poolsex male_vcf female_vcf output_file window min_reads range" << std::endl;
        exit(0);
    }

    std::string input_file_path = argv[1];
    std::ifstream input_file;
    input_file.open(input_file_path);

    std::string output_file_path = argv[2];
    std::ofstream output_file;
    output_file.open(output_file_path);

    uint32_t window = std::stoi(argv[3]);
    uint16_t min_reads_sex = std::stoi(argv[4]);
    float range = std::stof(argv[5]);

    if (not input_file.is_open()) {

        std::cout << "Error: cannot open input file." << std::endl;
        exit(0);
    }

    if (not output_file.is_open()) {

        std::cout << "Error: cannot open output file." << std::endl;
        exit(0);
    }

    std::cout << "Analyzing" << std::endl;
    analysis(input_file, output_file, window, min_reads_sex, range);

    return 0;
}
