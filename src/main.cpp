#include <iostream>
#include "variant_data.h"
#include "vcf_parsing.h"
#include "output.h"


int main(int argc, char *argv[]) {

    if (argc < 5) {
        std::cout << "Usage: poolsex male_vcf female_vcf output_file window" << std::endl;
        exit(0);
    }
    std::string male_file_path = argv[1];
    std::ifstream male_vcf_file;
    male_vcf_file.open(male_file_path);

    std::string female_file_path = argv[2];
    std::ifstream female_vcf_file;
    female_vcf_file.open(female_file_path);

    std::string output_file_path = argv[3];
    std::ofstream output_file;
    output_file.open(output_file_path);

//    uint16_t window = std::stoi(argv[4]);

    if (not male_vcf_file.is_open()) {

        std::cout << "Error: cannot open male VCF file." << std::endl;
        exit(0);
    }

    if (not female_vcf_file.is_open()) {

        std::cout << "Error: cannot open female VCF file." << std::endl;
        exit(0);
    }

    if (not output_file.is_open()) {

        std::cout << "Error: cannot open output file." << std::endl;
        exit(0);
    }

//    std::cout << "Getting contig lengths" << std::endl;
//    std::unordered_map<std::string, uint32_t> contig_lengths;
//    get_contig_lengths(male_vcf_file, contig_lengths);

//    std::cout << "Creating data structure" << std::endl;

//    variants data;
//    std::vector<Position> temp;
//    for (auto contig : contig_lengths) {
//        data[contig.first] = temp;
//        data[contig.first].resize(contig.second);
//    }

    std::cout << "Filling data structure" << std::endl;
    get_variant_data(male_vcf_file, female_vcf_file);

//    std::cout << "Generating output" << std::endl;
//    output_data(data, output_file);

    return 0;
}
