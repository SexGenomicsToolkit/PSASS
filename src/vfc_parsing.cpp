#include "vcf_parsing.h"



void get_contig_lengths(std::ifstream& vcf_file, std::unordered_map<std::string, uint32_t>& contig_lengths) {

    std::string line = "";
    std::vector<std::string> fields;
    std::string name = "";
    uint32_t size = 0;

    while(std::getline(vcf_file, line) and (line.substr(0, 2) == "##")) {
        if (line.substr(0, 8) == "##contig") {
            name = "";
            size = 0;
            fields = split(line, "=");
            name = split(fields[2], ",")[0];
            size = std::stoi(fields[3].substr(0, fields[3].size()-1));
            contig_lengths[name] = size;
        }
    }
}




void get_variant_data(std::ifstream& file, variants& data) {

    std::string line = "";
    std::string contig = "";
    uint64_t position = 0;
    size_t start, pos;

    while(std::getline(file, line)) {

        if (line.substr(0, 1) != "#" and line.size() > 1) {

            start = pos + 1;
            pos = line.substr(start).find("\t");
            position = std::stoi(line.substr(start, pos));

            std::cout << contig << " | " << position << std::endl;

//            temp = split(fields[9], ":");
//            coverage = std::stoi(temp[1]);
//            main_coverage = std::stoi(split(temp[2], ",")[0]);
//            data[contig][position].male_coverage += male_coverage;
//            if (male_coverage > 0 and main_coverage / male_coverage < 0.9) {
//                data[contig][position].male_heterozygote += 1;
//            }
        }
    }
}



//void get_variant_data(std::ifstream& file, variants& data, uint16_t window, std::unordered_map<std::string, uint32_t>& contig_lengths) {

//    variants data;
//    std::vector<Position> temp;
//    for (auto contig : contig_lengths) {
//        data[contig.first] = temp;
//        data[contig.first].resize(std::ceil(double(contig.second) / double(window)));
//    }

//    std::string line = "";
//    std::vector<std::string> fields, temp;
//    std::string contig = "";
//    uint64_t position = 0;
//    uint16_t male_coverage, main_coverage, leftover;

//    while(std::getline(file, line)) {
//        if (line.substr(0, 1) != "#") {
//            fields.resize(0);
//            fields = split(line, "\t");
//            contig = fields[0];
//            position = std::floor(stoi(fields[1]) / window);
//            temp = split(fields[9], ":");
//            male_coverage = std::stoi(temp[1]);
//            main_coverage = std::stoi(split(temp[2], ",")[0]);
//            data[contig][position].male_coverage += male_coverage;
//            if (male_coverage > 0 and main_coverage / male_coverage < 0.9) {
//                data[contig][position].male_heterozygote += 1;
//            }
//        }
//    }

//    for (auto& contig: data) {
//        for (uint i=0; i<contig.second.size()-1; ++i) {
//            contig.second[i].male_coverage /= window;
//        }
//        leftover = contig_lengths[contig.first] % window;
//        if (leftover != 0) {
//            contig.second[contig.second.size()-1].male_coverage /= leftover;
//        } else {
//            contig.second[contig.second.size()-1].male_coverage /= window;
//        }
//    }
//}

