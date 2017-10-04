#include <iostream>
#include "variant_data.h"
#include "vcf_parsing.h"
#include "output.h"


int main(int argc, char *argv[]) {

    std::string file_path = argv[1];
    std::ifstream vcf_file;
    vcf_file.open(file_path);

    std::string output_file_path = argv[2];
    std::ofstream output_file;
    output_file.open(output_file_path);

    uint16_t window = std::stoi(argv[3]);
    std::cout << window << std::endl;

    if (not vcf_file.is_open()) {

        std::cout << "Error: cannot open RADSeq VCF file." << std::endl;
        exit(0);
    }

    if (not output_file.is_open()) {

        std::cout << "Error: cannot open output file." << std::endl;
        exit(0);
    }

    std::cout << "Getting contig lengths" << std::endl;
    std::unordered_map<std::string, uint32_t> contig_lengths;
    get_contig_lengths(vcf_file, contig_lengths);

    std::cout << "Creating data structure" << std::endl;

    variants data;
    std::vector<Position> temp;
    for (auto contig : contig_lengths) {
        data[contig.first] = temp;
        data[contig.first].resize(std::ceil(double(contig.second) / double(window)));
    }

    std::cout << "Filling data structure" << std::endl;
    get_variant_data(vcf_file, data, window, contig_lengths);

    std::cout << "Generating output" << std::endl;
    output_data(data, output_file);

    return 0;
}
