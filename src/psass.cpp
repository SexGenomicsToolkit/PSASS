#include "psass.h"

// Psass class constructor
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

    this->output_handler = OutputHandler(&this->parameters, &this->input_data, this->male_pool, this->female_pool, this->male_index, this->female_index);
}



// Check whether current position is a sex-specific SNPs for each sex and update window_base_data.snps
void Psass::update_snps() {

    if (this->window_base_data.depth[0] > this->parameters.min_depth and this->window_base_data.depth[1] > this->parameters.min_depth) {

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
}



// Update window_base_data.coverage from each pool's data depth
void Psass::update_depth() {

    // Update data to push in window
    this->window_base_data.depth[0] = this->pair_data.pool1.depth;
    this->window_base_data.depth[1] = this->pair_data.pool2.depth;

    // Update total depth count to compute relative coverage later
    if (this->parameters.output_depth) {
        this->total_depth[0] += this->window_base_data.depth[0];
        this->total_depth[1] += this->window_base_data.depth[1];
    }
}



// Update sliding window data
void Psass::update_window() {

    // Add the current base data to the window (reset window if size bigger than window_size, which should not happen)
    (this->window.data.size() <= this->parameters.window_size) ? this->window.data.push_back(this->window_base_data) : this->window.data.resize(0);

    // If the window has size window_size, compute sum for each metric (first time with complete window in the current contig)
    if (this->window.data.size() <= this->parameters.window_size) {

        this->window.snps_in_window[0] += this->window_base_data.snps[0];
        this->window.snps_in_window[1] += this->window_base_data.snps[1];
        this->window.depth_in_window[0] += this->window_base_data.depth[0];
        this->window.depth_in_window[1] += this->window_base_data.depth[1];

    } else if (this->window.data.size() == this->parameters.window_size + 1) {  // Normal case (within contig) : substract front, add new value, remove front.

        this->window.snps_in_window[0] = this->window.snps_in_window[0] - this->window.data[0].snps[0] + this->window_base_data.snps[0];
        this->window.snps_in_window[1] = this->window.snps_in_window[1] - this->window.data[0].snps[1] + this->window_base_data.snps[1];
        this->window.depth_in_window[0] = this->window.depth_in_window[0] - this->window.data[0].depth[0] + this->window_base_data.depth[0];
        this->window.depth_in_window[1] = this->window.depth_in_window[1] - this->window.data[0].depth[1] + this->window_base_data.depth[1];
        this->window.data.pop_front();
    }
}



// Function called on a line from the input file (i.e. when meeting a '\n')
void Psass::process_line() {

    // Fill last pool2 base
    this->pair_data.pool2.nucleotides[5] = fast_stoi(this->input_data.temp.c_str());

    // Reset values
    this->pair_data.fst = 0;
    this->window_base_data.snps[0] = false;
    this->window_base_data.snps[1] = false;

    // Update data (depth per pool, fst, pi ...)
    this->pair_data.update();

    // Update depth data for window
    this->update_depth();

    // Update SNPs data for window
    this->update_snps();

    // Update window data
    this->update_window();

    ++this->total_bases;

    // Output SNPs positions
    if (parameters.output_snps_pos) {

        if (this->window_base_data.snps[this->male_index]) this->output_handler.output_snp_position("M");
        if (this->window_base_data.snps[this->female_index]) this->output_handler.output_snp_position("F");
    }

    // Output window information and update coverage
    if ((this->input_data.position - this->parameters.window_range) % this->parameters.output_resolution == 0 and this->input_data.position >= this->parameters.window_range) {

        if (parameters.output_snps_win) this->output_handler.output_snp_window(this->window.snps_in_window);

        if (parameters.output_depth) {

            this->depth_data[this->input_data.contig][this->input_data.position - this->parameters.window_range][0] = this->window.depth_in_window[this->male_index];
            this->depth_data[this->input_data.contig][this->input_data.position - this->parameters.window_range][1] = this->window.depth_in_window[this->female_index];
        }
    }

    // Change of contig
    if (this->input_data.contig != this->input_data.current_contig) {

        if (this->input_data.current_contig != "") {

            std::cout << "Finished analyzing contig :  " << this->input_data.current_contig << std::endl;

            if (parameters.output_snps_win) {
                this->window.snps_in_window[0] = 0;
                this->window.snps_in_window[1] = 0;
            }

            if (parameters.output_depth) {
                this->window.depth_in_window[0] = 0;
                this->window.depth_in_window[1] = 0;
            }

            this->window.data.resize(0);
        }
    }

    this->input_data.current_contig = this->input_data.contig;
    this->input_data.field = 0;
    this->input_data.temp = "";
}



// Function called on a field from the input file (i.e. when meeting a '\t')
void Psass::process_field() {

    switch (this->input_data.field) {

    case 0:
        this->input_data.contig = this->input_data.temp;
        break;

    case 1:
        this->input_data.position = fast_stoi(this->input_data.temp.c_str());
        break;

    case 2:
        break;

    case 3:
        this->pair_data.pool1.nucleotides[5] = fast_stoi(this->input_data.temp.c_str());
        break;

    default:
        break;
    }

    this->input_data.temp = "";
    this->input_data.subfield = 0;
    ++this->input_data.field;
}



// Function called on a subfield from the input file (i.e. when meeting a ':')
void Psass::process_subfield() {

    switch (this->input_data.field) {

    case 3:
        this->pair_data.pool1.nucleotides[this->input_data.subfield] = fast_stoi(this->input_data.temp.c_str());
        break;

    case 4:
        this->pair_data.pool2.nucleotides[this->input_data.subfield] = fast_stoi(this->input_data.temp.c_str());
        break;

    default:
        break;
    }
    this->input_data.temp = "";
    ++this->input_data.subfield;
}



// Read the input file and process each line
void Psass::run() {

    do {

        this->parameters.input_file.read(this->input_data.buff, this->input_data.buff_size);
        this->input_data.k = this->parameters.input_file.gcount();

        for (uint i=0; i<this->input_data.k; ++i) {

            switch (this->input_data.buff[i]) {

                case '\r':
                    break;

                case '\n':
                    this->process_line();
                    break;

                case '\t':
                    this->process_field();
                    break;

                case ':':
                    this->process_subfield();
                    break;

                default:
                    this->input_data.temp += this->input_data.buff[i];
                    break;

            }
        }

    } while (parameters.input_file);

    if (this->parameters.output_depth) this->output_handler.output_depth(this->depth_data, this->total_depth, this->total_bases);

}
