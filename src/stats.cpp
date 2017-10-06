#include "stats.h"

void calculate_fst(variants& data, std::ofstream& output_file, uint32_t window) {

    output_file << "Contig" << "\t" << "Position" << "\t" << "Fst" << "\n";

    float fst_window = 0;

    // Total reads
    float males_total = 0, females_total = 0;

    // Males / Females frequency for adenine, thymine, guanine, cytosine and indels
    float males_a_freq = 0, males_t_freq = 0, males_g_freq = 0, males_c_freq = 0, males_i_freq = 0;
    float females_a_freq = 0, females_t_freq = 0, females_g_freq = 0, females_c_freq = 0, females_i_freq = 0;

    // Average frequency for adenine, thymine, guanine, cytosine and indels
    float average_a_freq = 0, average_t_freq = 0, average_g_freq = 0, average_c_freq = 0, average_i_freq = 0;

    // Pi and Fst
    float males_pi = 0, females_pi = 0, total_pi = 0, within_pi = 0, fst = 0;

    for (auto contig : data) {

        for (uint i=0; i<contig.second.size(); ++i) {

            // Males and females totals
            males_total = float(contig.second[i].male_bases[0] + contig.second[i].male_bases[1] + contig.second[i].male_bases[2] +
                          contig.second[i].male_bases[3] + contig.second[i].male_bases[4]);
            females_total = float(contig.second[i].female_bases[0] + contig.second[i].female_bases[1] + contig.second[i].female_bases[2] +
                            contig.second[i].female_bases[3] + contig.second[i].female_bases[4]);

            if (males_total > 1 and females_total > 1) {

                // Male and females frequencies
                males_a_freq = float(contig.second[i].male_bases[0] / males_total);
                males_t_freq = float(contig.second[i].male_bases[1] / males_total);
                males_g_freq = float(contig.second[i].male_bases[2] / males_total);
                males_c_freq = float(contig.second[i].male_bases[3] / males_total);
                males_i_freq = float(contig.second[i].male_bases[4] / males_total);
                females_a_freq = float(contig.second[i].female_bases[0] / females_total);
                females_t_freq = float(contig.second[i].female_bases[1] / females_total);
                females_g_freq = float(contig.second[i].female_bases[2] / females_total);
                females_c_freq = float(contig.second[i].female_bases[3] / females_total);
                females_i_freq = float(contig.second[i].female_bases[4] / females_total);

                // Average frequencies
                average_a_freq = float((males_a_freq + females_a_freq) / 2);
                average_t_freq = float((males_t_freq + females_t_freq) / 2);
                average_g_freq = float((males_g_freq + females_g_freq) / 2);
                average_c_freq = float((males_c_freq + females_c_freq) / 2);
                average_i_freq = float((males_i_freq + females_i_freq) / 2);

                // Males, females, and total pi
                males_pi = float((1 - std::pow(males_a_freq, 2) - std::pow(males_t_freq, 2) - std::pow(males_g_freq, 2) - std::pow(males_c_freq, 2) - std::pow(males_i_freq, 2)) *
                           (males_total / (males_total - 1)));
                females_pi = float((1 - std::pow(females_a_freq, 2) - std::pow(females_t_freq, 2) - std::pow(females_g_freq, 2) - std::pow(females_c_freq, 2) - std::pow(females_i_freq, 2)) *
                           (females_total / (females_total - 1)));
                total_pi = float((1 - std::pow(average_a_freq, 2) - std::pow(average_t_freq, 2) - std::pow(average_g_freq, 2) - std::pow(average_c_freq, 2) - std::pow(average_i_freq, 2)) *
                           (std::min(males_total, females_total) / (std::min(males_total, females_total) - 1)));

                // Within pi
                within_pi = float((males_pi + females_pi) / 2);

                // Fst
                if (total_pi > 0) {
                    fst = float((total_pi - within_pi) / total_pi);
                } else {
                    fst = 0;
                }

                // Debug logs
//                if(contig.first == "LG1" and i > 3710 and i < 3720){
//                std::cout << contig.first << "\t" << i+1 << "\t" << males_total << "\t" << females_total << " | "
//                          << contig.second[i].male_bases[0] << "\t" << contig.second[i].male_bases[1] << "\t" << contig.second[i].male_bases[2] << "\t" << contig.second[i].male_bases[3] << "\t" << contig.second[i].male_bases[4] << " | "
//                          << contig.second[i].female_bases[0] << "\t" << contig.second[i].female_bases[1] << "\t" << contig.second[i].female_bases[2] << "\t" << contig.second[i].female_bases[3] << "\t" << contig.second[i].female_bases[4] << " | "
//                          << males_a_freq << "\t" << males_c_freq << "\t" << males_g_freq << "\t" << males_t_freq << "\t" << males_i_freq << "\t" << " | " << "\t"
//                          << females_a_freq << "\t" << females_c_freq << "\t" << females_g_freq << "\t" << females_t_freq << "\t" << females_i_freq << "\t" << " | " << "\t"
//                          << average_a_freq << "\t" << average_c_freq << "\t" << average_g_freq << "\t" << average_t_freq << "\t" << average_i_freq << "\t" << " | " << "\t"
//                          << males_pi << "\t" << females_pi << "\t" << total_pi << "\t" <<  within_pi << "\t" << fst << "\n";
//                }

            } else if (males_total > 1 or females_total > 1){

                fst = 1;

            } else {

                fst = 0;

            }

            // Output for a given window
            if (i % window == 0 and i > 0) {

                output_file << contig.first << "\t" << i / window - 1 << "\t" << float(fst_window / window) << "\n";
                fst_window = 0;

            } else if (i == contig.second.size() - 1) {

                output_file << contig.first << "\t" << uint(i / window) << "\t" << float(fst_window / (i % window)) << "\n";
                fst_window = 0;
            }

            fst_window += fst;

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
