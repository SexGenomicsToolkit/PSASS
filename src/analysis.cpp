#include "analysis.h"

void analysis(Parameters& parameters) {

    const uint window_range = parameters.window_size / 2;

    uint16_t pool1_bases[6], pool2_bases[6];
    bool snp_1 = false, snp_2 = false;

    // Total reads
    float pool1_total = 0, pool2_total = 0;

    // pool1 / pool2 frequency for adenine, thymine, guanine, cytosine and indels
    float pool1_a_freq = 0, pool1_t_freq = 0, pool1_g_freq = 0, pool1_c_freq = 0, pool1_i_freq = 0;
    float pool2_a_freq = 0, pool2_t_freq = 0, pool2_g_freq = 0, pool2_c_freq = 0, pool2_i_freq = 0;

    // Average frequency for adenine, thymine, guanine, cytosine and indels
    float average_a_freq = 0, average_t_freq = 0, average_g_freq = 0, average_c_freq = 0, average_i_freq = 0;

    // Pi and Fst
    float pool1_pi = 0, pool2_pi = 0, total_pi = 0, within_pi = 0, fst = 0;

    // Reading optimization parameters
    char buff[2048];
    uint k=0, field=0, subfield=0, position = 0;
    std::string contig = "", current_contig = "";
    std::string temp = "";

    // Sliding window
    uint fst_window = 0, snp_1_window = 0, snp_2_window = 0;
    std::deque<bool> fst_sliding_window, snp_1_sliding_window, snp_2_sliding_window;

    float coverage_1_window = 0, coverage_2_window = 0;
    std::deque<float> coverage_1_sliding_window, coverage_2_sliding_window;
    uint64_t total_coverage_1 = 0, total_coverage_2 = 0;
    table coverage;

    uint64_t total_bases = 0;

    do {

        parameters.input_file.read(buff, sizeof(buff));
        k = parameters.input_file.gcount();

        for (uint i=0; i<k; ++i) {

            switch (buff[i]) {

            case '\r':
                break;

            case '\n':
                // Fill last pool2 base
                pool2_bases[5] = fast_stoi(temp.c_str());

                // Reset values
                fst = 0;
                snp_1 = false;
                snp_2 = false;

                // pool1 and pool2 totals
                pool1_total = float(pool1_bases[0] + pool1_bases[1] + pool1_bases[2] + pool1_bases[3] + pool1_bases[5]);
                pool2_total = float(pool2_bases[0] + pool2_bases[1] + pool2_bases[2] + pool2_bases[3] + pool2_bases[5]);

                if (pool1_total > parameters.min_depth and pool2_total > parameters.min_depth) {

                    // pool1 and pool2 frequencies
                    pool1_a_freq = float(pool1_bases[0] / pool1_total);
                    pool1_t_freq = float(pool1_bases[1] / pool1_total);
                    pool1_g_freq = float(pool1_bases[2] / pool1_total);
                    pool1_c_freq = float(pool1_bases[3] / pool1_total);
                    pool1_i_freq = float(pool1_bases[5] / pool1_total);
                    pool2_a_freq = float(pool2_bases[0] / pool2_total);
                    pool2_t_freq = float(pool2_bases[1] / pool2_total);
                    pool2_g_freq = float(pool2_bases[2] / pool2_total);
                    pool2_c_freq = float(pool2_bases[3] / pool2_total);
                    pool2_i_freq = float(pool2_bases[5] / pool2_total);

                    // Average frequencies
                    average_a_freq = float((pool1_a_freq + pool2_a_freq) / 2);
                    average_t_freq = float((pool1_t_freq + pool2_t_freq) / 2);
                    average_g_freq = float((pool1_g_freq + pool2_g_freq) / 2);
                    average_c_freq = float((pool1_c_freq + pool2_c_freq) / 2);
                    average_i_freq = float((pool1_i_freq + pool2_i_freq) / 2);

                    // pool1, pool2, and total pi
                    pool1_pi = float((1 - pool1_a_freq * pool1_a_freq - pool1_t_freq * pool1_t_freq - pool1_g_freq * pool1_g_freq - pool1_c_freq * pool1_c_freq - pool1_i_freq * pool1_i_freq) *
                                     (pool1_total / (pool1_total - 1)));
                    pool2_pi = float((1 - pool2_a_freq * pool2_a_freq - pool2_t_freq * pool2_t_freq - pool2_g_freq * pool2_g_freq - pool2_c_freq * pool2_c_freq - pool2_i_freq * pool2_i_freq) *
                                       (pool2_total / (pool2_total - 1)));
                    total_pi = float((1 - average_a_freq * average_a_freq - average_t_freq * average_t_freq - average_g_freq * average_g_freq - average_c_freq * average_c_freq - average_i_freq * average_i_freq) *
                                     (std::min(pool1_total, pool2_total) / (std::min(pool1_total, pool2_total) - 1)));

                    // Within pi
                    within_pi = float((pool1_pi + pool2_pi) / 2);

                    // Fst
                    if (total_pi > 0) {
                        fst = float((total_pi - within_pi) / total_pi);
                    } else {
                        fst = 0;
                    }

                    // pool1 specific snps
                    if ((pool1_a_freq > 0.5 - parameters.snp_range and pool1_a_freq < 0.5 + parameters.snp_range and pool2_a_freq > 1 - parameters.fixed_range) or
                            (pool1_t_freq > 0.5 - parameters.snp_range and pool1_t_freq < 0.5 + parameters.snp_range and pool2_t_freq > 1 - parameters.fixed_range) or
                            (pool1_g_freq > 0.5 - parameters.snp_range and pool1_g_freq < 0.5 + parameters.snp_range and pool2_g_freq > 1 - parameters.fixed_range) or
                            (pool1_c_freq > 0.5 - parameters.snp_range and pool1_c_freq < 0.5 + parameters.snp_range and pool2_c_freq > 1 - parameters.fixed_range) or
                            (pool1_i_freq > 0.5 - parameters.snp_range and pool1_i_freq < 0.5 + parameters.snp_range and pool2_i_freq > 1 - parameters.fixed_range)) {

                        snp_1 = true;
                    }

                    // pool2 specific snps
                    if ((pool2_a_freq > 0.5 - parameters.snp_range and pool2_a_freq < 0.5 + parameters.snp_range and pool1_a_freq > 1 - parameters.fixed_range) or
                             (pool2_t_freq > 0.5 - parameters.snp_range and pool2_t_freq < 0.5 + parameters.snp_range and pool1_t_freq > 1 - parameters.fixed_range) or
                             (pool2_g_freq > 0.5 - parameters.snp_range and pool2_g_freq < 0.5 + parameters.snp_range and pool1_g_freq > 1 - parameters.fixed_range) or
                             (pool2_c_freq > 0.5 - parameters.snp_range and pool2_c_freq < 0.5 + parameters.snp_range and pool1_c_freq > 1 - parameters.fixed_range) or
                             (pool2_i_freq > 0.5 - parameters.snp_range and pool2_i_freq < 0.5 + parameters.snp_range and pool1_i_freq > 1 - parameters.fixed_range)) {

                         snp_2 = true;
                     }

                } // Could handle other cases, they could be interesting too

                // Refactor this part
                if (fst > parameters.min_fst) {
                    parameters.fst_threshold_output_file << contig << "\t" << position << "\t" << fst << "\n";
                }

                if (fst_sliding_window.size() <= parameters.window_size) {
                    fst_sliding_window.push_back(fst > parameters.min_fst);
                } else {
                    fst_sliding_window.resize(0);
                }
                if (fst_sliding_window.size() == parameters.window_size) {
                    fst_window = std::accumulate(fst_sliding_window.begin(), fst_sliding_window.end(), 0);
                } else if (fst_sliding_window.size() == parameters.window_size + 1) {
                    fst_window -= fst_sliding_window[0];
                    fst_window += (fst > parameters.min_fst);
                    fst_sliding_window.pop_front();
                }

                if (snp_1_sliding_window.size() <= parameters.window_size) {
                    snp_1_sliding_window.push_back(snp_1);
                } else {
                    snp_1_sliding_window.resize(0);
                }
                if (snp_1_sliding_window.size() == parameters.window_size) {
                    snp_1_window = 1.0*std::accumulate(snp_1_sliding_window.begin(), snp_1_sliding_window.end(), 0.0);
                } else if (snp_1_sliding_window.size() == parameters.window_size + 1) {
                    snp_1_window -= snp_1_sliding_window[0];
                    snp_1_window += snp_1;
                    snp_1_sliding_window.pop_front();
                }

                if (snp_2_sliding_window.size() <= parameters.window_size) {
                    snp_2_sliding_window.push_back(snp_2);
                } else {
                    snp_2_sliding_window.resize(0);
                }
                if (snp_2_sliding_window.size() == parameters.window_size) {
                    snp_2_window = 1.0*std::accumulate(snp_2_sliding_window.begin(), snp_2_sliding_window.end(), 0.0);
                } else if (snp_2_sliding_window.size() == parameters.window_size + 1) {
                    snp_2_window -= snp_2_sliding_window[0];
                    snp_2_window += snp_2;
                    snp_2_sliding_window.pop_front();
                }

                if (parameters.output_coverage) {

                    if (coverage_1_sliding_window.size() <= parameters.window_size) {
                        coverage_1_sliding_window.push_back(pool1_total);
                    } else {
                        coverage_1_sliding_window.resize(0);
                    }
                    if (coverage_1_sliding_window.size() == parameters.window_size) {
                        coverage_1_window = 1.0*std::accumulate(coverage_1_sliding_window.begin(), coverage_1_sliding_window.end(), 0.0);
                    } else if (coverage_1_sliding_window.size() == parameters.window_size + 1) {
                        coverage_1_window -= coverage_1_sliding_window[0];
                        coverage_1_window += pool1_total;
                        coverage_1_sliding_window.pop_front();
                    }

                    if (coverage_2_sliding_window.size() <= parameters.window_size) {
                        coverage_2_sliding_window.push_back(pool2_total);
                    } else {
                        coverage_2_sliding_window.resize(0);
                    }
                    if (coverage_2_sliding_window.size() == parameters.window_size) {
                        coverage_2_window = 1.0*std::accumulate(coverage_2_sliding_window.begin(), coverage_2_sliding_window.end(), 0.0);
                    } else if (coverage_2_sliding_window.size() == parameters.window_size + 1) {
                        coverage_2_window -= coverage_2_sliding_window[0];
                        coverage_2_window += pool2_total;
                        coverage_2_sliding_window.pop_front();
                    }
                    total_coverage_1 += uint(pool1_total);
                    total_coverage_2 += uint(pool2_total);
                }

                if (parameters.output_snps_pos) {

                    if (snp_1 and snp_2) std::cout << "SNP in both sexes !" << std::endl;

                    if (snp_1) {

                        if (parameters.male_pool == 1) {

                            parameters.snps_pos_output_file << contig << "\t" << position << "\t" << "M" << "\n";

                        } else {

                            parameters.snps_pos_output_file << contig << "\t" << position << "\t" << "F" << "\n";

                        }

                    } else if (snp_2) {

                        if (parameters.male_pool == 1) {

                            parameters.snps_pos_output_file << contig << "\t" << position << "\t" << "F" << "\n";

                        } else {

                            parameters.snps_pos_output_file << contig << "\t" << position << "\t" << "M" << "\n";

                        }

                    }
                }

                ++total_bases;

                if (contig != current_contig and current_contig != "") {

                    std::cout << "Finished analyzing contig :  " << current_contig << std::endl;

                    if ((position - window_range) % parameters.output_resolution == 0 and position > parameters.window_size) {
                        parameters.fst_window_output_file << contig << "\t" << position - window_range << "\t" << fst_window << "\n";
                        parameters.snps_output_file << contig << "\t" << position - window_range << "\t";
                        if (parameters.male_pool == 1) {
                            parameters.snps_output_file << snp_1_window << "\t" << snp_2_window << "\n";
                            if (parameters.output_coverage) coverage[contig][position] = std::pair<float, float>(coverage_1_window, coverage_2_window);
                        } else {
                            parameters.snps_output_file << snp_2_window << "\t" << snp_1_window << "\n";
                            if (parameters.output_coverage) coverage[contig][position] = std::pair<float, float>(coverage_2_window, coverage_1_window);
                        }
                    }

                    fst_window = 0;
                    snp_1_window = 0;
                    snp_2_window = 0;
                    fst_sliding_window.resize(0);
                    snp_1_sliding_window.resize(0);
                    snp_2_sliding_window.resize(0);

                    if (parameters.output_coverage) {
                        coverage_2_window = 0;
                        coverage_1_window = 0;
                        coverage_1_sliding_window.resize(0);
                        coverage_2_sliding_window.resize(0);
                    }

                } else {

                    if ((position - window_range) % parameters.output_resolution == 0 and position > parameters.window_size) {
                        parameters.fst_window_output_file << contig << "\t" << position - window_range << "\t" << fst_window << "\n";
                        parameters.snps_output_file << contig << "\t" << position - window_range << "\t";
                        if (parameters.male_pool == 1) {
                            parameters.snps_output_file << snp_1_window << "\t" << snp_2_window << "\n";
                            if (parameters.output_coverage) coverage[contig][position] = std::pair<float, float>(coverage_1_window, coverage_2_window);
                        } else {
                            parameters.snps_output_file << snp_2_window << "\t" << snp_1_window << "\n";
                            if (parameters.output_coverage) coverage[contig][position] = std::pair<float, float>(coverage_2_window, coverage_1_window);
                        }
                    }
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
                    position = fast_stoi(temp.c_str());

                case 2:
                    break;

                case 3:
                    pool1_bases[5] = fast_stoi(temp.c_str());
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
                    pool1_bases[subfield] = fast_stoi(temp.c_str());
                    break;

                case 4:
                    pool2_bases[subfield] = fast_stoi(temp.c_str());
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

    } while (parameters.input_file);

    if (fst > parameters.min_fst) {
        parameters.fst_threshold_output_file << contig << "\t" << position << "\t" << fst << "\n";
    }

    if ((position - window_range) % parameters.output_resolution == 0 and position > parameters.window_size) {
        parameters.fst_window_output_file << contig << "\t" << position - window_range << "\t" << fst_window << "\n";
        parameters.snps_output_file << contig << "\t" << position - window_range << "\t";
        if (parameters.male_pool == 1) {
            parameters.snps_output_file << snp_1_window << "\t" << snp_2_window << "\n";
        } else {
            parameters.snps_output_file << snp_2_window << "\t" << snp_1_window << "\n";
        }
    }

    if (parameters.output_coverage) {

        float average_coverage_m = 0;
        float average_coverage_f = 0;
        if (parameters.male_pool == 1) {
            average_coverage_m = float(total_coverage_1) / float(total_bases);
            average_coverage_f = float(total_coverage_2) / float(total_bases);
        } else {
            average_coverage_m = float(total_coverage_2) / float(total_bases);
            average_coverage_f = float(total_coverage_1) / float(total_bases);
        }

        for (auto const& contig : coverage) {
            for (auto const& position: contig.second) {
                    parameters.coverage_output_file << contig.first << "\t" << position.first << "\t" << std::fixed << std::setprecision(2) <<
                                                       float((position.second.first / parameters.window_size)/ average_coverage_m) <<
                                                       "\t" << float((position.second.second / parameters.window_size) / average_coverage_f) <<
                                                       "\t" << float(position.second.first / parameters.window_size) <<
                                                       "\t" << float(position.second.second / parameters.window_size) << "\n";
            }
        }
    }
}
