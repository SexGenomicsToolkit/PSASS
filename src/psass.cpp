#include "psass.h"

// Psass class constructor
Psass::Psass(int argc, char *argv[]) {

    this->t_begin = std::chrono::steady_clock::now();


    ArgParser cmd_options(argc, argv);
    cmd_options.set_parameters(this->parameters);

    this->logs = Logs(this->parameters.output_prefix);
    this->logs.write("PSASS started.");
    cmd_options.output_parameters(logs.file);

    this->output_handler = OutputHandler(&this->parameters, &this->input_data, &this->pair_data, &this->depth_data, &this->gff_data.genes, &this->logs);

    std::cout << "Preprocessing data ..." << std::endl;
    if (this->parameters.output_genes) {

        this->gff_data.read_gff_file(this->parameters.gff_file, this->logs);
        this->parameters.output_genes = true;

    }
}



// Update window_base_data.nucleotides
void Psass::update_fst_parts() {

    // Fst computation for the window is implemented using the formula described in Karlsson et al 2007
    // https://www.nature.com/articles/ng.2007.10   https://doi.org/10.1038/ng.2007.10

    uint8_t n_alleles[2];
    uint32_t allele_p1 = 0, allele_p2 = 0;
    float a1 = 0.0, a2 = 0.0, n1 = 0.0, n2 = 0.0;
    float h1 = 0.0, h2 = 0.0, N = 0.0, D = 0.0;

    if (this->pair_data.pool1.depth > this->parameters.min_depth and this->pair_data.pool2.depth > this->parameters.min_depth) {  // Only compute FST when each pool has depth higher than min depth threshold

        n_alleles[0] = 0;
        n_alleles[1] = 0;
        allele_p1 = 0;
        allele_p2 = 0;
        a1 = 0;
        a2 = 0;
        n1 = 0;
        n2 = 0;
        h1 = 0;
        h2 = 0;
        N = 0;
        D = 0;

        for (uint i=0; i<5; ++i) {  // Count the number of alleles in each pool (do not include indels)

            if (this->pair_data.pool1.nucleotides[i] > 0) {

                ++n_alleles[0];
                allele_p1 = this->pair_data.pool1.nucleotides[i];
                allele_p2 = this->pair_data.pool2.nucleotides[i];

            }

            if (this->pair_data.pool2.nucleotides[i] > 0) {

                ++n_alleles[1];
            }
        }

        if (n_alleles[0] == 2 and n_alleles[1] == 2) {  // Only compute FST for biallelic positions

            a1 = float(allele_p1);
            a2 = float(allele_p2);
            n1 = float(this->pair_data.pool1.depth);
            n2 = float(this->pair_data.pool2.depth);

            h1 = a1 * (n1 - a1) / (n1 * (n1 - 1));
            h2 = a2 * (n2 - a2) / (n2 * (n2 - 1));
            N = (a1 / n1 - a2 / n2) * (a1 / n1 - a2 / n2) - h1 / n1 - h2 / n2;
            D = N + h1 + h2;

        }
    }

    this->window_base_data.fst_parts[0] = N;
    this->window_base_data.fst_parts[1] = D;

}



// Check whether current position is a sex-specific SNPs for each sex and update window_base_data.snps
void Psass::update_snps() {

    this->window_base_data.snps[0] = false;
    this->window_base_data.snps[1] = false;

    if (this->pair_data.pool1.depth > this->parameters.min_depth and this->pair_data.pool2.depth > this->parameters.min_depth) {

        for (auto i=0; i<6; ++i) {

            if (this->pair_data.pool1.frequencies[i] > this->parameters.min_het and
                this->pair_data.pool1.frequencies[i] < this->parameters.max_het and
                this->pair_data.pool2.frequencies[i] > this->parameters.min_hom) {

                if (this->parameters.group_snps) {

                    if (not this->consecutive_snps[0]) this->window_base_data.snps[0] = true;
                    else this->window_base_data.snps[0] = false;
                    this->consecutive_snps[0] = true;

                } else {

                    this->window_base_data.snps[0] = true;

                }

            }

            if (this->pair_data.pool2.frequencies[i] > this->parameters.min_het and
                this->pair_data.pool2.frequencies[i] < this->parameters.max_het and
                this->pair_data.pool1.frequencies[i] > this->parameters.min_hom) {

                if (this->parameters.group_snps) {

                    if (not this->consecutive_snps[1]) this->window_base_data.snps[1] = true;
                    else this->window_base_data.snps[1] = false;
                    this->consecutive_snps[1] = true;

                } else {

                    this->window_base_data.snps[1] = true;

                }

            }
        }

        if (not this->window_base_data.snps[0] and not this->window_base_data.snps[1]) {

            this->consecutive_snps[0] = false;
            this->consecutive_snps[1] = false;

        }
    }
}



// Update window_base_data.coverage from each pool's data depth
void Psass::update_depth() {

    // Update data to push in window
    if (this->parameters.output_depth) {

        this->window_base_data.depth[0] = this->pair_data.pool1.depth;
        this->window_base_data.depth[1] = this->pair_data.pool2.depth;

    }

    // Update total depth count to compute relative coverage later
    if (this->parameters.output_depth or this->parameters.output_genes) {

        this->total_depth[0] += this->pair_data.pool1.depth;
        this->total_depth[1] += this->pair_data.pool2.depth;

    }
}



// Update sliding window data
void Psass::update_window(bool end) {

    // Add the current base data to the window (reset window if size bigger than window_size, which should not happen)
    (this->window.data.size() <= this->parameters.window_size and not end) ? this->window.data.push_back(this->window_base_data) : this->window.data.resize(0);

    // If the window is smaller than window_size, only add to total (beginning of the contig)
    if (this->window.data.size() <= this->parameters.window_size) {

        if (this->parameters.output_snps_win) {

            this->window.snps_in_window[0] += this->window_base_data.snps[0];
            this->window.snps_in_window[1] += this->window_base_data.snps[1];

        }

        if (this->parameters.output_depth) {

            this->window.depth_in_window[0] += this->window_base_data.depth[0];
            this->window.depth_in_window[1] += this->window_base_data.depth[1];

        }

        if (this->parameters.output_fst_win) {

            this->window.fst_parts[0] += this->window_base_data.fst_parts[0];
            this->window.fst_parts[1] += this->window_base_data.fst_parts[1];

        }


    } else if (this->window.data.size() == this->parameters.window_size + 1) {  // Normal case (within contig) : substract front, add new value, remove front.

        if (this->parameters.output_fst_win) {

            this->window.fst_parts[0] = this->window.fst_parts[0]
                                        - this->window.data[0].fst_parts[0]
                                        + this->window_base_data.fst_parts[0];

            this->window.fst_parts[1] = this->window.fst_parts[1]
                                        - this->window.data[0].fst_parts[1]
                                        + this->window_base_data.fst_parts[1];

        }

        if (this->parameters.output_snps_win) {

            this->window.snps_in_window[0] = this->window.snps_in_window[0]
                                                            - this->window.data[0].snps[0]
                                                            + this->window_base_data.snps[0];

            this->window.snps_in_window[1] = this->window.snps_in_window[1]
                                                              - this->window.data[0].snps[1]
                                                              + this->window_base_data.snps[1];

        }

        if (this->parameters.output_depth) {

            this->window.depth_in_window[0] = this->window.depth_in_window[0]
                                                             - this->window.data[0].depth[0]
                                                             + this->pair_data.pool1.depth;

            this->window.depth_in_window[1] = this->window.depth_in_window[1]
                                                               - this->window.data[0].depth[1]
                                                               + this->pair_data.pool2.depth;

        }

        this->window.data.pop_front();

    }
}



// Update genes data
void Psass::update_genes() {

    std::string gene = "", contig = "";
    bool coding = false;

    if (this->gff_data.contig.find(this->input_data.position) != this->gff_data.contig.end()) {

        gene = this->gff_data.contig[this->input_data.position].first;
        coding = this->gff_data.contig[this->input_data.position].second;

        // Coding or non-coding depth and snps
        this->gff_data.genes[gene].depth[2 * 0 + coding] += this->pair_data.pool1.depth;
        this->gff_data.genes[gene].depth[2 * 1 + coding] += this->pair_data.pool2.depth;
        this->gff_data.genes[gene].snps[2 * 0 + coding] += this->window_base_data.snps[0];
        this->gff_data.genes[gene].snps[2 * 1 + coding] += this->window_base_data.snps[1];

        // Gene-level depth and snps
        this->gff_data.genes[gene].depth[4 + 0] += this->pair_data.pool1.depth;
        this->gff_data.genes[gene].depth[4 + 1] += this->pair_data.pool2.depth;
        this->gff_data.genes[gene].snps[4 + 0] += this->window_base_data.snps[0];
        this->gff_data.genes[gene].snps[4 + 1] += this->window_base_data.snps[1];

    }
}



// Handles output of sliding window
void Psass::output_window_step() {

    if (this->parameters.output_fst_win) this->output_handler.output_fst_window(this->window.fst_parts);

    if (this->parameters.output_snps_win) this->output_handler.output_snp_window(this->window.snps_in_window);

    if (this->parameters.output_depth) {

        this->depth_data[this->input_data.contig][this->input_data.position - this->parameters.window_range][0] = this->window.depth_in_window[0];
        this->depth_data[this->input_data.contig][this->input_data.position - this->parameters.window_range][1] = this->window.depth_in_window[1];
        this->depth_data[this->input_data.contig][this->input_data.position - this->parameters.window_range][2] = float(this->window.data.size());

    }
}



// End of contig needs special processing. Progressively remove the beginning of the window until last position.
void Psass::process_contig_end() {

    uint last_spot = uint(this->input_data.last_position / this->parameters.output_resolution) * this->parameters.output_resolution + this->parameters.window_range;
    uint first_spot = last_spot - this->parameters.window_range;
    auto tmp = this->window.data;

    uint32_t tmp_position = this->input_data.position;
    std::string tmp_contig = this->input_data.contig;

    this->input_data.contig = this->input_data.current_contig;

    for (auto i = first_spot; i <= last_spot; i += this->parameters.output_resolution) {

        this->window.fst_parts[0] = 0;
        this->window.fst_parts[1] = 0;
        this->window.snps_in_window[0] = 0;
        this->window.snps_in_window[1] = 0;
        this->window.depth_in_window[0] = 0;
        this->window.depth_in_window[1] = 0;

        this->input_data.position = i;

        for (uint j = 0; j < this->parameters.output_resolution / 2; ++j) tmp.pop_front();

        for (auto base: tmp) {

            this->window.snps_in_window[0] += base.snps[0];
            this->window.snps_in_window[1] += base.snps[1];
            this->window.depth_in_window[0] += base.depth[0];
            this->window.depth_in_window[1] += base.depth[1];
            this->window.fst_parts[0] += base.fst_parts[0];
            this->window.fst_parts[1] += base.fst_parts[1];
        }

        this->window.data = tmp;

        this->output_window_step();

    }

    this->input_data.position = tmp_position;
    this->input_data.contig = tmp_contig;
}



// Function called on a line from the input file (i.e. when meeting a '\n')
void Psass::process_line() {

    // Fill last pool2 base
    this->pair_data.pool2.nucleotides[5] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));

    // Reset values
    this->pair_data.fst = 0;
    this->window_base_data.snps[0] = false;
    this->window_base_data.snps[1] = false;

    // Update / reset values if change of contig
    if (this->input_data.contig != this->input_data.current_contig) {

        if (this->input_data.current_contig != "") {

            if(this->window.data.size() >= this->parameters.window_size) this->process_contig_end();

            this->logs.write("Processing of contig <" + this->input_data.current_contig + "> ended without errors.");
            this->logs.write("Processing of contig <" + this->input_data.contig + "> started.");

            if (parameters.output_snps_win) {

                this->window.snps_in_window[0] = 0;
                this->window.snps_in_window[1] = 0;

            }

            if (parameters.output_depth) {

                this->window.depth_in_window[0] = 0;
                this->window.depth_in_window[1] = 0;

            }

            if (parameters.output_fst_win) {

                this->window.fst_parts[0] = 0;
                this->window.fst_parts[1] = 0;

            }

            this->window.data.resize(0);

        } else {

            this->logs.write("Processing of contig <" + this->input_data.contig + "> started.");

        }

        if (parameters.output_genes) this->gff_data.new_contig(this->input_data, this->logs);

        if (parameters.group_snps) {

            this->consecutive_snps[0] = false;
            this->consecutive_snps[1] = false;

        }
    }

    // Reset line parsing values
    this->input_data.current_contig = this->input_data.contig;
    this->input_data.last_position = this->input_data.position;
    this->input_data.field = 0;
    this->input_data.temp = "";

    // Update data (depth per pool, fst, pi ...)
    this->pair_data.update(this->parameters.output_fst_pos);

    // Update nucleotides data for window
    if (this->parameters.output_fst_win) this->update_fst_parts();

    // Update depth data for window
    if (this->parameters.output_depth or this->parameters.output_genes) this->update_depth();

    // Update SNPs data for window
    if (this->parameters.output_snps_win or this->parameters.output_snps_pos or this->parameters.output_genes) this->update_snps();

    // Update window data
    if (this->parameters.output_fst_win or this->parameters.output_snps_win or this->parameters.output_depth) this->update_window();

    // Update genes data
    if (this->parameters.output_genes) this->update_genes();

    ++this->total_bases;

    // Output Fst positions
    if (parameters.output_fst_pos) {

        if (this->pair_data.fst > this->parameters.min_fst) this->output_handler.output_fst_position(this->pair_data.fst);
    }

    // Output SNPs positions
    if (parameters.output_snps_pos) {

        if (this->window_base_data.snps[0]) this->output_handler.output_snp_position(parameters.pool1_id);
        if (this->window_base_data.snps[1]) this->output_handler.output_snp_position(parameters.pool2_id);

    }

    // Output window information and update coverage
    if ((this->input_data.position == 1 or (this->input_data.position - this->parameters.window_range) % this->parameters.output_resolution == 0)
        and this->input_data.position >= this->parameters.window_range) {

        this->output_window_step();

    }
}



// Function called on a field from the input file (i.e. when meeting a '\t')
void Psass::process_popoolation_field() {

    switch (this->input_data.field) {

    case 0:  // Scaffold field
        this->input_data.contig = this->input_data.temp;
        break;

    case 1:  // Position field
        this->input_data.position = static_cast<uint32_t>(fast_stoi(this->input_data.temp.c_str()));
        break;

    case 2:  // Reference base
        break;

    case 3:  // Frequencies in first pool (last subfield is indels)
        this->pair_data.pool1.nucleotides[5] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
        break;

    default:
        break;

    }

    this->input_data.temp = "";
    this->input_data.subfield = 0;
    ++this->input_data.field;
}



// Function called on a subfield from the input file (i.e. when meeting a ':')
void Psass::process_popoolation_subfield() {

    switch (this->input_data.field) {

    case 3:
        this->pair_data.pool1.nucleotides[this->input_data.subfield] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
        break;

    case 4:
        this->pair_data.pool2.nucleotides[this->input_data.subfield] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
        break;

    default:
        break;
    }
    this->input_data.temp = "";
    ++this->input_data.subfield;
}



// Function called on a field from the input file (i.e. when meeting a '\t')
void Psass::process_psass_field() {

    switch (this->input_data.field) {

        case 0:
            this->input_data.contig = this->input_data.temp;
            break;

        case 1:
            this->input_data.position = static_cast<uint32_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 2:
            break;

        case 3:
            this->pair_data.pool1.nucleotides[0] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 4:
            this->pair_data.pool1.nucleotides[1] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 5:
            this->pair_data.pool1.nucleotides[2] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 6:
            this->pair_data.pool1.nucleotides[3] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 7:
            this->pair_data.pool1.nucleotides[4] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 8:
            this->pair_data.pool1.nucleotides[5] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 9:
            this->pair_data.pool2.nucleotides[0] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 10:
            this->pair_data.pool2.nucleotides[1] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 11:
            this->pair_data.pool2.nucleotides[2] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 12:
            this->pair_data.pool2.nucleotides[3] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 13:
            this->pair_data.pool2.nucleotides[4] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        case 14:
            this->pair_data.pool2.nucleotides[5] = static_cast<uint16_t>(fast_stoi(this->input_data.temp.c_str()));
            break;

        default:
            break;

    }

    this->input_data.temp = "";
    ++this->input_data.field;
}



// Read the input file and process each line
void Psass::run() {

    this->logs.write("Processing of <" + this->parameters.input_file_path + "> started.");
    std::cerr << "PSASS started." << std::endl;

    do {

        this->parameters.input_file.read(this->input_data.buff, this->input_data.buff_size);
        this->input_data.k = this->parameters.input_file.gcount();

        if (this->parameters.popoolation_format) {

            for (uint i=0; i<this->input_data.k; ++i) {

                switch (this->input_data.buff[i]) {

                    case ':':
                        this->process_popoolation_subfield();
                        break;

                    case '\t':
                        this->process_popoolation_field();
                        break;

                    case '\n':
                        this->process_line();
                        break;

                    default:
                        this->input_data.temp += this->input_data.buff[i];
                        break;

                }
            }

        } else {

            for (uint i=0; i<this->input_data.k; ++i) {

                switch (this->input_data.buff[i]) {

                    case '\t':
                        this->process_psass_field();
                        break;

                    case '\n':
                        this->process_line();
                        break;

                    default:
                        this->input_data.temp += this->input_data.buff[i];
                        break;

                }
            }

        }

    } while (parameters.input_file);

    this->input_data.contig = "";
    this->process_line();

    this->logs.write("Processing of <" + this->parameters.input_file_path + "> ended without errors.");
    this->logs.write("Processed <" + std::to_string(this->total_bases) + "> lines.");  // One base per line

    this->average_depth[0] = float(this->total_depth[0]) / float(this->total_bases);
    this->average_depth[1] = float(this->total_depth[1]) / float(this->total_bases);

    if (this->parameters.output_depth) this->output_handler.output_depth(this->average_depth);
    if (this->parameters.output_genes) this->output_handler.output_genes(this->average_depth);

    std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    long seconds = std::chrono::duration_cast<std::chrono::seconds>(t_end - t_begin).count();
    long minutes = seconds / 60;
    long hours = minutes / 60;
    this->logs.write("Total runtime : " + std::to_string(hours) + "h " + std::to_string(minutes%60) + "m " + std::to_string(seconds%60) + "s.");
    this->logs.write("PSASS ended without errors.");
    std::cerr << "PSASS ended successfully." << std::endl;

}
