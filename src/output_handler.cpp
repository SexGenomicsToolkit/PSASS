#include "output_handler.h"

OutputHandler::OutputHandler(Parameters* parameters, InputData* input_data, PoolBaseData* male_pool, PoolBaseData* female_pool, bool male_index, bool female_index) {

    // Pointers to data structures from PSASS
    this->input_data = input_data;
    this->male_pool = male_pool;
    this->female_pool = female_pool;
    this->male_index = male_index;
    this->female_index = female_index;

    this->parameters = parameters;

    // Create base output file path
    if (this->parameters->output_prefix != "") this->parameters->output_prefix += "_";

    // Open output file objects
    if (this->parameters->output_fst_pos) this->fst_position_output_file.open(this->parameters->output_prefix);
    if (this->parameters->output_fst_win) this->fst_window_output_file.open(this->parameters->output_prefix);
    if (this->parameters->output_snps_pos) this->snps_position_output_file.open(this->parameters->output_prefix);
    if (this->parameters->output_snps_win) this->snps_window_output_file.open(this->parameters->output_prefix);
    if (this->parameters->output_coverage) this->coverage_output_file.open(this->parameters->output_prefix);
    if (this->parameters->output_genes) this->genes_output_file.open(this->parameters->output_prefix);
}

// Write SNP and nucleotide information if current base is a sex-specific SNP
void OutputFile::open(const std::string& prefix) {

    this->path = prefix + this->suffix + ".tsv";
    this->file.open(this->path);
    if (not this->file.is_open()) {
        std::cerr << "Error: cannot open output file (" <<this->path << ")." << std::endl;
        exit(1);
    }
    this->file << this->header;
}



// Write SNP and nucleotide information if current base is a sex-specific SNP
void OutputHandler::output_snp_position(std::string sex) {

    this->snps_position_output_file.file << std::fixed << std::setprecision(2)
                                         << this->input_data->contig << "\t" << this->input_data->position << "\t" << sex << "\t"
                                         << this->male_pool->output_frequencies() << "\t"
                                         << this->female_pool->output_frequencies() << "\n";
}



// Write SNP and nucleotide information if current base is a sex-specific SNP
void OutputHandler::output_snp_window(uint16_t snps_total[2]) {

    this->snps_window_output_file.file << this->input_data->contig << "\t" << this->input_data->position - this->parameters->window_range << "\t"
                                       << snps_total[this->male_index] << "\t"
                                       << snps_total[this->female_index] << "\n";
}

