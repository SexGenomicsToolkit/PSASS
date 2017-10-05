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

//struct Position {
//    uint16_t male_coverage = 0;
//    uint16_t female_coverage = 0;
//    bool male_heterozygote = 0;
//    bool female_heterozygote = 0;
//    bool male_specific_heterozygote = 0;
//    bool female_specific_heterozygote = 0;
//    uint16_t male_bases[4] {0, 0, 0, 0};
//    uint16_t female_bases[4] {0, 0, 0, 0};
//    double fst = 0;
//    double fisher_p = 0;
//    double cochran_p = 0;
//};


void get_variant_data(std::ifstream& male_file, std::ifstream& female_file, std::output& output_file, uint16_t window) {

    std::string male_line = "", female_line = "";
    std::string male_contig = "", female_contig = "";
    uint64_t male_position = 0, female_position = 0;
    size_t male_start, male_pos, female_start, female_pos;
    uint16_t males[5] {0, 0, 0, 0, 0}, females[5] {0, 0, 0, 0, 0};
    bool male_catchup = false, female_catchup = false;

    while(male_file and female_file) {

        if (!male_catchup) std::getline(male_file, male_line);
        if (!female_catchup) std::getline(female_file, female_line);

        if (male_line.substr(0, 1) != "#" and male_line.size() > 1) {

            male_pos = male_line.find("\t");
            male_contig = male_line.substr(0, male_pos);
            male_start = male_pos + 1;
            male_pos = male_line.substr(male_start).find("\t");
            male_position = std::stoi(male_line.substr(male_start, male_pos));

        }

        if (female_line.substr(0, 1) != "#" and female_line.size() > 1) {

            female_pos = female_line.find("\t");
            female_contig = female_line.substr(0, female_pos);
            female_start = female_pos + 1;
            female_pos = female_line.substr(female_start).find("\t");
            female_position = std::stoi(female_line.substr(female_start, female_pos));
        }

        if (male_contig == female_contig) {

            if (male_position == female_position) {

                output_file << male_contig << "\t" << male_position << "\t" << 0 << std::endl;

            } else if (male_position < female_position) {

                male_catchup = true;
                output_file << male_contig << "\t" << male_position << "\t" << 1 << std::endl;
            }

        if (male_contig != "" and female_contig != "") std::cout << male_contig << " | " << male_position << "  |  " <<
                                                                    female_contig << " | " << female_position << std::endl;

//            start = pos + 1;
//            pos = line.substr(start).find("\t");
//            position = std::stoi(line.substr(start, pos));

//            std::cout << contig << " | " << position << std::endl;

//            temp = split(fields[9], ":");
//            coverage = std::stoi(temp[1]);
//            main_coverage = std::stoi(split(temp[2], ",")[0]);
//            data[contig][position].male_coverage += male_coverage;
//            if (male_coverage > 0 and main_coverage / male_coverage < 0.9) {
//                data[contig][position].male_heterozygote += 1;
//            }
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

