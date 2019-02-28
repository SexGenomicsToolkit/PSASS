#include "output_handler.h"

OutputHandler::OutputHandler(Parameters* parameters, InputData* input_data, PoolBaseData* male_pool, PoolBaseData* female_pool, bool male_index, bool female_index,
                             std::map<std::string, std::map<uint, float[3]>>* depth, std::unordered_map<std::string, Gene>* genes, Logs* logs) {

    // Pointers to data structures from PSASS
    this->input_data = input_data;
    this->male_pool = male_pool;
    this->female_pool = female_pool;
    this->male_index = male_index;
    this->female_index = female_index;
    this->depth = depth;
    this->genes = genes;
    this->logs = logs;

    this->parameters = parameters;

    // Create base output file path
    if (this->parameters->output_prefix != "") this->parameters->output_prefix += "_";

    // Open output file objects
    if (this->parameters->output_fst_pos) this->fst_position_output_file.open(this->parameters->output_prefix, this->logs);
    if (this->parameters->output_fst_win) this->fst_window_output_file.open(this->parameters->output_prefix, this->logs);
    if (this->parameters->output_snps_pos) this->snps_position_output_file.open(this->parameters->output_prefix, this->logs);
    if (this->parameters->output_snps_win) this->snps_window_output_file.open(this->parameters->output_prefix, this->logs);
    if (this->parameters->output_depth) this->depth_output_file.open(this->parameters->output_prefix, this->logs);
    if (this->parameters->output_genes) this->genes_output_file.open(this->parameters->output_prefix, this->logs);
}



// Write SNP and nucleotide information if current base is a sex-specific SNP
void OutputFile::open(const std::string& prefix, Logs* logs) {

    this->path = prefix + this->suffix + ".tsv";
    this->file.open(this->path);

    if (not this->file.is_open()) {

        std::cerr << "Error: cannot open output file <" << this->path << ">." << std::endl;
        logs->write("Error: cannot open output file <" + this->path + ">.");
        exit(1);

    }

    logs->write("Created output file <" + this->path + ">.");
    this->file << this->header;
}



// Write Fst information if current base has fst higher than specified threshold
void OutputHandler::output_fst_position(float fst) {

    this->fst_position_output_file.file << std::fixed << std::setprecision(4)
                                        << this->input_data->contig << "\t" << this->input_data->position << "\t" << fst <<  "\n";
}



// Write Fst information for the current windows
void OutputHandler::output_fst_window(float fst_parts[2]) {

    float fst = 0;

    (fst_parts[1] > 0) ? fst = std::max(float(0.0), fst_parts[0] / fst_parts[1]) : fst = 0.0;

    this->fst_window_output_file.file << this->input_data->contig << "\t" << this->input_data->position - this->parameters->window_range << "\t"
                                      << std::fixed << std::setprecision(4)
                                      << fst << "\n";
}



// Write SNP and nucleotide information if current base is a sex-specific SNP
void OutputHandler::output_snp_position(std::string sex) {

    this->snps_position_output_file.file << this->input_data->contig << "\t" << this->input_data->position << "\t" << sex << "\t"
                                         << this->male_pool->output_frequencies() << "\t"
                                         << this->female_pool->output_frequencies() << "\n";
}



// Write SNP and nucleotide information for the current window
void OutputHandler::output_snp_window(uint32_t snps_total[2]) {

    this->snps_window_output_file.file << this->input_data->contig << "\t" << this->input_data->position - this->parameters->window_range << "\t"
                                       << snps_total[this->male_index] << "\t"
                                       << snps_total[this->female_index] << "\n";
}


// Write depth information at the end of the analysis
void OutputHandler::output_depth(float* average_depth) {

    logs->write("Depth data output started.");

    float window_size = 0;

    for (auto const& contig : *this->depth) {

        for (auto const& position: contig.second) {

            window_size = position.second[2];

            this->depth_output_file.file << contig.first << "\t" << position.first << "\t"
                                         << float(position.second[this->male_index] / window_size) << "\t"
                                         << float(position.second[this->female_index] / window_size) << "\t"
                                         << std::fixed << std::setprecision(2)
                                         << float((position.second[this->male_index] / window_size)/ average_depth[this->male_index]) << "\t"
                                         << float((position.second[this->female_index] / window_size)/ average_depth[this->female_index]) << "\n";

        }
    }

    logs->write("Depth data output ended without errors.");
}


// Write genes information at the end of the analysis
void OutputHandler::output_genes(float* average_depth) {

    logs->write("Genes data output started.");

    float depth_correction_males = (average_depth[this->male_index] + average_depth[this->female_index]) / 2 / average_depth[this->male_index];
    float depth_correction_females = (average_depth[this->male_index] + average_depth[this->female_index]) / 2 / average_depth[this->female_index];

    uint gene_length = 0, male_depth = 0, female_depth = 0;

    for (auto gene: *this->genes) {

        gene_length = uint(std::stoi(gene.second.end) -  std::stoi(gene.second.start));
        male_depth = (gene.second.depth[4 + this->male_index]) / gene_length;
        female_depth = (gene.second.depth[4 + this->female_index]) / gene_length;
        gene.second.noncoding_length = gene_length - gene.second.coding_length;

        if (gene.second.coding_length == 0) {

            gene.second.depth[2 * this->male_index] = 0;
            gene.second.depth[2 * this->female_index] = 0;

        } else {

            (gene.second.noncoding_length > 0) ? gene.second.depth[2 * this->male_index] /= gene.second.noncoding_length : gene.second.depth[2 * this->male_index] = 0;
            (gene.second.noncoding_length > 0) ? gene.second.depth[2 * this->female_index] /= gene.second.noncoding_length: gene.second.depth[2 * this->female_index] = 0;

        }

        if (gene.second.noncoding_length == 0) {

            gene.second.depth[2 * this->male_index + 1] = 0;
            gene.second.depth[2 * this->female_index + 1] = 0;

        } else {

            (gene.second.coding_length > 0) ? gene.second.depth[2 * this->male_index + 1] /= gene.second.coding_length : gene.second.depth[2 * this->male_index + 1] = 0;
            (gene.second.coding_length > 0) ? gene.second.depth[2 * this->female_index + 1] /= gene.second.coding_length: gene.second.depth[2 * this->female_index + 1] = 0;

        }

        this->genes_output_file.file << gene.second.contig << "\t" << gene.second.start << "\t" << gene.second.end << "\t"
                                     <<gene.second.id << "\t" << gene.second.name << "\t" << gene.second.product << "\t"
                                     << male_depth << "\t" << int(male_depth * depth_correction_males) << "\t"
                                     << gene.second.depth[2 * this->male_index + 1] << "\t" << int(gene.second.depth[2 * this->male_index + 1] * depth_correction_males) << "\t"
                                     << gene.second.depth[2 * this->male_index] << "\t" << int(gene.second.depth[2 * this->male_index] * depth_correction_males) << "\t"
                                     << female_depth << "\t" << int(female_depth * depth_correction_females) << "\t"
                                     << gene.second.depth[2 * this->female_index + 1] << "\t" << int(gene.second.depth[2 * this->female_index + 1] * depth_correction_females) << "\t"
                                     << gene.second.depth[2 * this->female_index] << "\t" << int(gene.second.depth[2 * this->female_index] * depth_correction_females) << "\t"
                                     << gene.second.snps[4 + this->male_index] << "\t" << gene.second.snps[2 * this->male_index + 1] << "\t" << gene.second.snps[2 * this->male_index] << "\t"
                                     << gene.second.snps[4 + this->female_index] << "\t" << gene.second.snps[2 * this->female_index + 1] << "\t" << gene.second.snps[2 * this->female_index] << "\n";
    }

    logs->write("Genes data output ended without errors.");

}
