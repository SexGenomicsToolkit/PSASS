#include "analysis.h"
#include <numeric>

uint analysis(Parameters& parameters) {

    // GFF data structures
    std::unordered_map<std::string, std::vector<std::vector<std::string>>> gff_data;
    std::unordered_map<std::string, Gene> genes;
    std::unordered_map<uint, std::pair<std::string, bool>> regions;

    if (parameters.output_genes) {

        write_log("\n", parameters.log_file, false, false);
        write_log("Reading GFF file...", parameters.log_file, true, true);

        read_gff_file(parameters.gff_file, gff_data, genes);

        write_log("Finished reading GFF file : ", parameters.log_file, true, false);
        write_log(genes.size(), parameters.log_file, false, false);
        write_log(" genes found.", parameters.log_file, false, true);
    }

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
    uint k=0, field=0, subfield=0, position = 0, n_lines = 0;
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

    // Gene filtering
    std::string gene;
    std::vector<std::string> infos;
    bool coding;


    write_log("Processing input sync file ...", parameters.log_file, true, true);

    do {

        parameters.input_file.read(buff, sizeof(buff));
        k = parameters.input_file.gcount();

        for (uint i=0; i<k; ++i) {

            switch (buff[i]) {

            case '\r':
                break;

            case '\n':

                // Update log file
                ++n_lines;
                if (n_lines % 5000000 == 0) {
                    write_log("Processed ", parameters.log_file, true, false);
                    write_log(n_lines / 1000000, parameters.log_file, false, false);
                    write_log(" M. lines.", parameters.log_file, false, true);
                }

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
                    if (parameters.output_fst_pos or parameters.output_fst_win) {
                        if (total_pi > 0) {
                            fst = float((total_pi - within_pi) / total_pi);
                        } else {
                            fst = 0;
                        }
                    }

                    // Snps
                    if (parameters.output_snps_pos or parameters.output_snps_win or parameters.output_genes) {

                        // pool1 specific snps
                        if ((pool1_a_freq > 0.5 - parameters.range_het and pool1_a_freq < 0.5 + parameters.range_het and pool2_a_freq > 1 - parameters.range_hom) or
                                (pool1_t_freq > 0.5 - parameters.range_het and pool1_t_freq < 0.5 + parameters.range_het and pool2_t_freq > 1 - parameters.range_hom) or
                                (pool1_g_freq > 0.5 - parameters.range_het and pool1_g_freq < 0.5 + parameters.range_het and pool2_g_freq > 1 - parameters.range_hom) or
                                (pool1_c_freq > 0.5 - parameters.range_het and pool1_c_freq < 0.5 + parameters.range_het and pool2_c_freq > 1 - parameters.range_hom) or
                                (pool1_i_freq > 0.5 - parameters.range_het and pool1_i_freq < 0.5 + parameters.range_het and pool2_i_freq > 1 - parameters.range_hom)) {

                            snp_1 = true;
                        }

                        // pool2 specific snps
                        if ((pool2_a_freq > 0.5 - parameters.range_het and pool2_a_freq < 0.5 + parameters.range_het and pool1_a_freq > 1 - parameters.range_hom) or
                                 (pool2_t_freq > 0.5 - parameters.range_het and pool2_t_freq < 0.5 + parameters.range_het and pool1_t_freq > 1 - parameters.range_hom) or
                                 (pool2_g_freq > 0.5 - parameters.range_het and pool2_g_freq < 0.5 + parameters.range_het and pool1_g_freq > 1 - parameters.range_hom) or
                                 (pool2_c_freq > 0.5 - parameters.range_het and pool2_c_freq < 0.5 + parameters.range_het and pool1_c_freq > 1 - parameters.range_hom) or
                                 (pool2_i_freq > 0.5 - parameters.range_het and pool2_i_freq < 0.5 + parameters.range_het and pool1_i_freq > 1 - parameters.range_hom)) {

                             snp_2 = true;
                         }
                    }

                } // Could handle other cases, they could be interesting too

                // FST positions
                if (parameters.output_fst_pos and fst > parameters.min_fst) {
                    parameters.fst_pos_output_file << contig << "\t" << position << "\t" << fst << "\n";
                }

                // FST sliding window
                if (parameters.output_fst_win) {
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
                }

                // SNPs positions
                if (parameters.output_snps_pos) {

                    if (snp_1 and snp_2) std::cout << "SNP in both sexes !" << std::endl;

                    if (snp_1) {

                        if (parameters.male_pool == 1) {

                            parameters.snps_pos_output_file << std::fixed << std::setprecision(2)
                                                            << contig << "\t" << position << "\t" << "M" << "\t"
                                                            << pool1_a_freq << "\t" << pool1_t_freq << "\t"
                                                            << pool1_g_freq << "\t" << pool1_c_freq << "\t"
                                                            << pool1_i_freq << "\t"
                                                            << pool2_a_freq << "\t" << pool2_t_freq << "\t"
                                                            << pool2_g_freq << "\t" << pool2_c_freq << "\t"
                                                            << pool2_i_freq << "\n";

                        } else {

                            parameters.snps_pos_output_file << std::fixed << std::setprecision(2)
                                                            << contig << "\t" << position << "\t" << "F" << "\t"
                                                            << pool2_a_freq << "\t" << pool2_t_freq << "\t"
                                                            << pool2_g_freq << "\t" << pool2_c_freq << "\t"
                                                            << pool2_i_freq << "\t"
                                                            << pool1_a_freq << "\t" << pool1_t_freq << "\t"
                                                            << pool1_g_freq << "\t" << pool1_c_freq << "\t"
                                                            << pool1_i_freq << "\n";

                        }

                    } else if (snp_2) {

                        if (parameters.male_pool == 1) {

                            parameters.snps_pos_output_file << std::fixed << std::setprecision(2)
                                                            << contig << "\t" << position << "\t" << "F" << "\t"
                                                            << pool1_a_freq << "\t" << pool1_t_freq << "\t"
                                                            << pool1_g_freq << "\t" << pool1_c_freq << "\t"
                                                            << pool1_i_freq << "\t"
                                                            << pool2_a_freq << "\t" << pool2_t_freq << "\t"
                                                            << pool2_g_freq << "\t" << pool2_c_freq << "\t"
                                                            << pool2_i_freq << "\n";

                        } else {

                            parameters.snps_pos_output_file << std::fixed << std::setprecision(2)
                                                            << contig << "\t" << position << "\t" << "M" << "\t"
                                                            << pool2_a_freq << "\t" << pool2_t_freq << "\t"
                                                            << pool2_g_freq << "\t" << pool2_c_freq << "\t"
                                                            << pool2_i_freq << "\t"
                                                            << pool1_a_freq << "\t" << pool1_t_freq << "\t"
                                                            << pool1_g_freq << "\t" << pool1_c_freq << "\t"
                                                            << pool1_i_freq << "\n";

                        }
                    }
                }

                // SNPs window
                if (parameters.output_snps_win) {

                    // SNP 1 sliding window
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

                    // SNP 2 sliding window
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
                }

                // Coverage
                if (parameters.output_coverage or parameters.output_genes) {

                    if (coverage_1_sliding_window.size() <= parameters.window_size) {
                        coverage_1_sliding_window.push_back(pool1_total);
                    } else {
                        coverage_1_sliding_window.resize(0);
                    }

                    if (coverage_1_sliding_window.size() == parameters.window_size) {
                        coverage_1_window = 1.0 * std::accumulate(coverage_1_sliding_window.begin(), coverage_1_sliding_window.end(), 0.0);
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

                // Genes
                if (parameters.output_genes) {

                    if (regions.find(position) != regions.end()) {

                        gene = regions[position].first;
                        coding = regions[position].second;

                        if (coding) {
                            if (parameters.male_pool == 1) {
                                genes[gene].coverage[0] += pool1_total;
                                genes[gene].coverage[2] += pool2_total;
                                genes[gene].snps[0] += snp_1;
                                genes[gene].snps[2] += snp_2;
                            } else {
                                genes[gene].coverage[0] += pool2_total;
                                genes[gene].coverage[2] += pool1_total;
                                genes[gene].snps[0] += snp_2;
                                genes[gene].snps[2] += snp_1;
                            }
                        } else {
                            if (parameters.male_pool == 1) {
                                genes[gene].coverage[1] += pool1_total;
                                genes[gene].coverage[3] += pool2_total;
                                genes[gene].snps[1] += snp_1;
                                genes[gene].snps[3] += snp_2;
                            } else {
                                genes[gene].coverage[1] += pool2_total;
                                genes[gene].coverage[3] += pool1_total;
                                genes[gene].snps[1] += snp_2;
                                genes[gene].snps[3] += snp_1;
                            }
                        }
                    }
                }

                ++total_bases;

                // Output window results (fst, snps, coverage)
                if ((position - window_range) % parameters.output_resolution == 0 and position > window_range) {
                    if (parameters.output_fst_win) parameters.fst_win_output_file << contig << "\t" << position - window_range << "\t" << fst_window << "\n";
                    if (parameters.output_snps_win) parameters.snps_win_output_file << contig << "\t" << position - window_range << "\t";
                    if (parameters.male_pool == 1) {
                        if (parameters.output_snps_win) parameters.snps_win_output_file << snp_1_window << "\t" << snp_2_window << "\n";
                        if (parameters.output_coverage) coverage[contig][position] = std::pair<float, float>(coverage_1_window, coverage_2_window);
                    } else {
                        if (parameters.output_snps_win) parameters.snps_win_output_file << snp_2_window << "\t" << snp_1_window << "\n";
                        if (parameters.output_coverage) coverage[contig][position] = std::pair<float, float>(coverage_2_window, coverage_1_window);
                    }
                }

                // Change of contig
                if (contig != current_contig) {

                    regions.clear();

                    for (auto line: gff_data[contig]) {

                        infos = split(line[8], ";");
                        for (auto i: infos) {
                            if (i.substr(0, 5) == "gene=") gene = split(i, "=")[1];
                        }

                        if (line[2] == "gene") {
                            for (int i=std::stoi(line[3]); i < std::stoi(line[4]) + 1; ++i) regions[i] = std::pair<std::string, bool>(gene, false);
                        } else {
                            for (int i=std::stoi(line[3]); i < std::stoi(line[4]) + 1; ++i) regions[i] = std::pair<std::string, bool>(gene, true);
                        }
                    }

                    if (current_contig != "") {

                        std::cout << "Finished analyzing contig :  " << current_contig << std::endl;

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

    // Output coverage results
    write_log("Generating coverage output file...", parameters.log_file, true, true);
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
                    parameters.coverage_output_file << contig.first << "\t" << position.first - window_range << "\t" << std::fixed << std::setprecision(2) <<
                                                       float((position.second.first / parameters.window_size)/ average_coverage_m) <<
                                                       "\t" << float((position.second.second / parameters.window_size) / average_coverage_f) <<
                                                       "\t" << float(position.second.first / parameters.window_size) <<
                                                       "\t" << float(position.second.second / parameters.window_size) << "\n";
            }
        }
    }

    // Output genes results
    if (parameters.output_genes) {

        uint gene_length, male_coverage, female_coverage;
        write_log("Generating genes output file...", parameters.log_file, true, true);

        for (auto gene: genes) {

            gene_length = std::stoi(gene.second.end) -  std::stoi(gene.second.start);
            male_coverage = (gene.second.coverage[0] + gene.second.coverage[1]) / gene_length;
            female_coverage = (gene.second.coverage[2] + gene.second.coverage[3]) / gene_length;
            gene.second.noncoding_length = gene_length -  gene.second.coding_length;

            if (gene.second.coding_length == 0) {
                gene.second.coverage[0] = 0;
                gene.second.coverage[2] = 0;
            } else {
                gene.second.coverage[0] /= gene.second.coding_length;
                gene.second.coverage[2] /= gene.second.coding_length;
            }

            if (gene.second.noncoding_length == 0) {
                gene.second.coverage[1] = 0;
                gene.second.coverage[3] = 0;
            } else {
                gene.second.coverage[1] /= gene.second.noncoding_length;
                gene.second.coverage[3] /= gene.second.noncoding_length;
            }

            parameters.genes_output_file << gene.second.contig << "\t" << gene.second.start << "\t" << gene.second.end << "\t" << gene.second.name << "\t" << gene.second.product << "\t"
                                         << male_coverage << "\t" << gene.second.coverage[0] << "\t" << gene.second.coverage[1] << "\t"
                                         << female_coverage << "\t" << gene.second.coverage[2] << "\t" << gene.second.coverage[3] << "\t"
                                         << gene.second.snps[0] + gene.second.snps[1] << "\t" << gene.second.snps[0] << "\t" << gene.second.snps[1] << "\t"
                                         << gene.second.snps[2] + gene.second.snps[3] << "\t" << gene.second.snps[2] << "\t" << gene.second.snps[3] << "\n";
        }
    }


    return n_lines;
}
