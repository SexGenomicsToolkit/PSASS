#include "psass.h"

Psass::Psass(int argc, char *argv[]) {

    ArgParser cmd_options(argc, argv);
    cmd_options.set_parameters(this->parameters);
    cmd_options.print_parameters();

    if (this->parameters.male_pool == 1) {
        this->male_pool = &this->pair_data.pool1;
        this->female_pool = &this->pair_data.pool2;
    } else {
        this->male_pool = &this->pair_data.pool2;
        this->female_pool = &this->pair_data.pool1;
    }

    this->male_index = (this->parameters.male_pool == 2);  // 0 if male pool is first, 1 otherwise
    this->female_index = (this->parameters.male_pool == 1);  // 0 if male pool is second, 1 otherwise

}


void Psass::update_snps() {

    this->window_base_data.snps[0] = false;
    this->window_base_data.snps[1] = false;

    for (auto i=0; i<6; ++i) {

        if (this->pair_data.pool1.frequencies[i] > this->parameters.min_het and
            this->pair_data.pool1.frequencies[i] < this->parameters.max_het and
            this->pair_data.pool2.frequencies[i] > this->parameters.min_hom) {

            this->window_base_data.snps[0] = true;
        }

        if (this->pair_data.pool2.frequencies[i] > this->parameters.min_het and
            this->pair_data.pool2.frequencies[i] < this->parameters.max_het and
            this->pair_data.pool1.frequencies[i] > this->parameters.min_hom) {

            this->window_base_data.snps[1] = true;
        }
    }
}


void Psass::update_coverage() {

    this->window_base_data.coverage[0] = this->pair_data.pool1.depth;
    this->window_base_data.coverage[1] = this->pair_data.pool2.depth;
}


void Psass::update_window() {

    // Add the current base to the window (reset window if size bigger than window_size, which should not happen)
    (this->window.data.size() <= this->parameters.window_size) ? this->window.data.push_back(this->window_base_data) : this->window.data.resize(0);

    // If the window has size window_size, compute sum for each metric (first time with complete window in the current contig)
    if (this->window.data.size() == this->parameters.window_size) {

        for (auto base: this->window.data) {
            this->window.snps_total[0] += base.snps[0];
            this->window.snps_total[1] += base.snps[1];
            this->window.coverage_total[0] += base.coverage[0];
            this->window.coverage_total[1] += base.coverage[1];
        }

    } else if (this->window.data.size() == this->parameters.window_size + 1) {  // Normal case (within contig) : substract front, add new value, remove front.

        this->window.snps_total[0] = this->window.snps_total[0] - this->window.data[0].snps[0] + this->window_base_data.snps[0];
        this->window.snps_total[1] = this->window.snps_total[1] - this->window.data[0].snps[1] + this->window_base_data.snps[1];
        this->window.coverage_total[0] = this->window.snps_total[0] - this->window.data[0].coverage[0] + this->window_base_data.coverage[0];
        this->window.coverage_total[1] = this->window.snps_total[1] - this->window.data[0].coverage[1] + this->window_base_data.coverage[1];
        this->window.data.pop_front();
    }
}


void Psass::output_snp(std::string sex, std::string& contig, uint position) {

    this->parameters.snps_pos_output_file << std::fixed << std::setprecision(2)
                                         << contig << "\t" << position << "\t" << sex << "\t"
                                         << this->male_pool->output_frequencies() << "\t"
                                         << this->female_pool->output_frequencies() << "\n";
}


void Psass::run() {

    char buff[2048];
    long k = 0;
    uint field = 0, subfield = 0, position = 0, n_lines = 0;
    std::string contig = "", current_contig = "";
    std::string temp = "";

    do {
        this->parameters.input_file.read(buff, sizeof(buff));
        k = this->parameters.input_file.gcount();

        for (uint i=0; i<k; ++i) {

            switch (buff[i]) {

            case '\r':
                break;

            case '\n':

                // Fill last pool2 base
                this->pair_data.pool2.nucleotides[5] = fast_stoi(temp.c_str());

                // Reset values
                this->pair_data.fst = 0;
                this->window_base_data.snps[0] = false;
                this->window_base_data.snps[1] = false;

                this->pair_data.update();
                this->update_coverage();

                if (this->window_base_data.coverage[0] > this->parameters.min_depth and this->window_base_data.coverage[1] > this->parameters.min_depth) {

                    this->update_snps();

                } // Could handle other cases, they could be interesting too

                // SNPs positions
                if (parameters.output_snps_pos) {

                    if (this->window_base_data.snps[this->male_index]) this->output_snp("M", contig, position);
                    if (this->window_base_data.snps[this->female_index]) this->output_snp("F", contig, position);

                }

                this->update_window();

                ++this->total_bases;

                // Output window information and update coverage
                if ((position - this->parameters.window_range) % this->parameters.output_resolution == 0 and position > this->parameters.window_range) {

                    if (parameters.output_snps_win) {
                        this->parameters.snps_win_output_file << contig << "\t" << position - this->parameters.window_range << "\t"
                                                              << this->window.snps_total[this->male_index] << "\t"
                                                              << this->window.snps_total[this->female_index] << "\n";
                    }

                    if (parameters.output_coverage) {
                        this->coverage[contig][position][0] = this->window.coverage_total[this->male_index];
                        this->coverage[contig][position][1] = this->window.coverage_total[this->female_index];
                    }

                }

                // Change of contig
                if (contig != current_contig) {

                    if (current_contig != "") {

                        std::cout << "Finished analyzing contig :  " << current_contig << std::endl;

                        if (parameters.output_snps_win) {
                            this->window.snps_total[0] = 0;
                            this->window.snps_total[1] = 0;
                        }

                        if (parameters.output_coverage) {
                            this->window.coverage_total[0] = 0;
                            this->window.coverage_total[1] = 0;
                        }

                        this->window.data.resize(0);
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
                    this->pair_data.pool1.nucleotides[5] = fast_stoi(temp.c_str());
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
                    this->pair_data.pool1.nucleotides[subfield] = fast_stoi(temp.c_str());
                    break;

                case 4:
                    this->pair_data.pool2.nucleotides[subfield] = fast_stoi(temp.c_str());
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
}
