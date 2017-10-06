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




void get_variant_data(std::ifstream& file, variants& data, bool male) {

    std::string line = "", contig = "", base = "";
    std::vector<std::string> fields, coverages, alleles;
    uint64_t position = 0;
    uint16_t coverage = 0;

    while(std::getline(file, line)) {

        if (line.substr(0, 1) != "#" and line.size() > 1) {

            fields = split(line, "\t");
            contig = fields[0];
            position = std::stoi(fields[1]) - 1;

            // Coverage for each allele
            coverages = split(split(fields[9], ":")[2], ",");

            // Reference allele
            base = fields[3];
            coverage = std::stoi(coverages[0]);

            if (base.size() == 1) {
                if (male) {
                    switch(base[0]){
                        case 'A':
                            data[contig][position].male_bases[0] = coverage;
                            break;
                        case 'T':
                            data[contig][position].male_bases[1] = coverage;
                            break;
                        case 'G':
                            data[contig][position].male_bases[2] = coverage;
                            break;
                        case 'C':
                            data[contig][position].male_bases[3] = coverage;
                            break;
                        default:
                            break;
                    }
                } else {
                    switch(base[0]){
                        case 'A':
                            data[contig][position].female_bases[0] = coverage;
                            break;
                        case 'T':
                            data[contig][position].female_bases[1] = coverage;
                            break;
                        case 'G':
                            data[contig][position].female_bases[2] = coverage;
                            break;
                        case 'C':
                            data[contig][position].female_bases[3] = coverage;
                            break;
                        default:
                            break;
                    }
                }
            }

            // Other alleles
            alleles = split(fields[4], ",");

            for (uint i=0; i < alleles.size() - 1; ++i) {

                base = alleles[i];
                coverage = std::stoi(coverages[i+1]);

                if (base.size() == 1) {
                    if (male) {
                        switch(base[0]){
                            case 'A':
                                data[contig][position].male_bases[0] = coverage;
                                break;
                            case 'T':
                                data[contig][position].male_bases[1] = coverage;
                                break;
                            case 'G':
                                data[contig][position].male_bases[2] = coverage;
                                break;
                            case 'C':
                                data[contig][position].male_bases[3] = coverage;
                                break;
                            default:
                                break;
                        }
                    } else {
                        switch(base[0]){
                            case 'A':
                                data[contig][position].female_bases[0] = coverage;
                                break;
                            case 'T':
                                data[contig][position].female_bases[1] = coverage;
                                break;
                            case 'G':
                                data[contig][position].female_bases[2] = coverage;
                                break;
                            case 'C':
                                data[contig][position].female_bases[3] = coverage;
                                break;
                            default:
                                break;
                        }
                    }
                } else {
                    if (male) data[contig][position].male_bases[4] = coverage;
                    else data[contig][position].female_bases[4] = coverage;
                }
            }
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

