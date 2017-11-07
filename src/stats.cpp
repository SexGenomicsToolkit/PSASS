#include "stats.h"

void analysis(std::ifstream& input_file, std::ofstream& fst_output_file, std::ofstream& snps_output_file, uint16_t min_reads_sex, float range, float min_fst) {

    fst_output_file << "Contig" << "\t" << "Position" << "\t" << "Fst" << "\n";
    snps_output_file << "Contig" << "\t" << "Position" << "\t" << "Sex" << "\n";

    uint16_t male_bases[6], female_bases[6];
    bool snp_m = false, snp_f = false;

    // Total reads
    float males_total = 0, females_total = 0;

    // Males / Females frequency for adenine, thymine, guanine, cytosine and indels
    float males_a_freq = 0, males_t_freq = 0, males_g_freq = 0, males_c_freq = 0, males_i_freq = 0;
    float females_a_freq = 0, females_t_freq = 0, females_g_freq = 0, females_c_freq = 0, females_i_freq = 0;

    // Average frequency for adenine, thymine, guanine, cytosine and indels
    float average_a_freq = 0, average_t_freq = 0, average_g_freq = 0, average_c_freq = 0, average_i_freq = 0;

    // Pi and Fst
    float males_pi = 0, females_pi = 0, total_pi = 0, within_pi = 0, fst = 0;

    // Reading optimization parameters
    char buff[2048];
    uint k=0, field=0, subfield=0, position = 0;
    std::string contig = "", current_contig = "";
    std::string temp = "";

    do {

        input_file.read(buff, sizeof(buff));
        k = input_file.gcount();

        for (uint i=0; i<k; ++i) {

            switch (buff[i]) {

            case '\r':
                break;

            case '\n':
                fst = 0;
                snp_f = false;
                snp_m = false;

                female_bases[5] = std::stoi(temp);

                // Males and females totals
                males_total = float(male_bases[0] + male_bases[1] + male_bases[2] + male_bases[3] + male_bases[5]);
                females_total = float(female_bases[0] + female_bases[1] + female_bases[2] + female_bases[3] + female_bases[5]);

                if (males_total > min_reads_sex and females_total > min_reads_sex) {

                    // Male and females frequencies
                    males_a_freq = float(male_bases[0] / males_total);
                    males_t_freq = float(male_bases[1] / males_total);
                    males_g_freq = float(male_bases[2] / males_total);
                    males_c_freq = float(male_bases[3] / males_total);
                    males_i_freq = float(male_bases[5] / males_total);
                    females_a_freq = float(female_bases[0] / females_total);
                    females_t_freq = float(female_bases[1] / females_total);
                    females_g_freq = float(female_bases[2] / females_total);
                    females_c_freq = float(female_bases[3] / females_total);
                    females_i_freq = float(female_bases[5] / females_total);

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

                    if ((males_a_freq > 0.5 - range and males_a_freq < 0.5 + range and females_a_freq > 0.9) or
                            (males_t_freq > 0.5 - range and males_t_freq < 0.5 + range and females_t_freq > 0.9) or
                            (males_g_freq > 0.5 - range and males_g_freq < 0.5 + range and females_g_freq > 0.9) or
                            (males_c_freq > 0.5 - range and males_c_freq < 0.5 + range and females_c_freq > 0.9) or
                            (males_i_freq > 0.5 - range and males_i_freq < 0.5 + range and females_i_freq > 0.9)) {

                        snp_m = true;
                    }

                    if ((females_a_freq > 0.5 - range and females_a_freq < 0.5 + range and males_a_freq > 0.9) or
                             (females_t_freq > 0.5 - range and females_t_freq < 0.5 + range and males_t_freq > 0.9) or
                             (females_g_freq > 0.5 - range and females_g_freq < 0.5 + range and males_g_freq > 0.9) or
                             (females_c_freq > 0.5 - range and females_c_freq < 0.5 + range and males_c_freq > 0.9) or
                             (females_i_freq > 0.5 - range and females_i_freq < 0.5 + range and males_i_freq > 0.9)) {

                         snp_f = true;
                     }


                }

                if (fst >= min_fst) {

                    fst_output_file << contig << "\t" << position << "\t" << fst << "\n";

                }

                if (snp_m and snp_f) {

                    snps_output_file << contig << "\t" << position << "\t" << "B" << "\n";

                } else if (snp_m) {

                    snps_output_file << contig << "\t" << position << "\t" << "M" << "\n";

                } else if (snp_f) {

                    snps_output_file << contig << "\t" << position << "\t" << "F" << "\n";

                }

                current_contig = contig;
                field = 0;
                temp = "";
                break;

            case '\t':

                switch (field) {

                case 0:
                    contig = temp;
                    break;

                case 1:
                    position = std::stoi(temp);

                case 2:
                    break;

                case 3:
                    male_bases[5] = std::stoi(temp);
                    break;

                default:
                    break;
                }

                temp = "";
                subfield = 0;
                ++field;
                break;

            case ':':

                switch (field) {

                case 3:
                    male_bases[subfield] = std::stoi(temp);
                    break;

                case 4:
                    female_bases[subfield] = std::stoi(temp);
                    break;

                default:
                    break;
                }
                temp = "";
                ++subfield;
                break;

            default:
                temp += buff[i];
                break;
            }

        }

    } while (input_file);

    if (fst > min_fst) {

        fst_output_file << contig << "\t" << position << "\t" << fst << "\n";

    }

    if (snp_m and snp_f) {

        snps_output_file << contig << "\t" << position << "\t" << "B" << "\n";

    } else if (snp_m) {

        snps_output_file << contig << "\t" << position << "\t" << "M" << "\n";

    } else if (snp_f) {

        snps_output_file << contig << "\t" << position << "\t" << "F" << "\n";

    }

}