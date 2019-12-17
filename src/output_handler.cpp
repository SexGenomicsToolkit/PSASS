#include "output_handler.h"

OutputHandler::OutputHandler(Parameters* parameters, InputData* input_data, PairBaseData* pair_data, std::map<std::string, std::map<uint, float[3]>>* depth, std::unordered_map<std::string, Gene>* genes) {

    // Pointers to data structures from PSASS
    this->input_data = input_data;
    this->pair_data = pair_data;
    this->depth = depth;
    this->genes = genes;

    this->parameters = parameters;

    // Create base output file path
    if (this->parameters->output_prefix != "") this->parameters->output_prefix += "_";

    // Open output file objects
    if (this->parameters->output_fst_pos) {
        this->fst_position_output_file.open(this->parameters->output_prefix);
        this->fst_position_output_file.file << "Contig" << "\t" << "Position" << "\t" << "Fst" << "\n";
    }

    if (this->parameters->output_fst_win) {
        this->fst_window_output_file.open(this->parameters->output_prefix);
        this->fst_window_output_file.file << "Contig" << "\t" << "Position" << "\t" << "Fst" << "\n";
    }

    if (this->parameters->output_snps_pos) {
        this->snps_position_output_file.open(this->parameters->output_prefix);
        this->snps_position_output_file.file << "Contig" << "\t" << "Position" << "\t" << "Pool" << "\t" <<
                                             this->parameters->pool1_id << "_A" << "\t" << this->parameters->pool1_id << "_T" << "\t" <<
                                             this->parameters->pool1_id << "_C" << "\t" << this->parameters->pool1_id << "_G" << "\t" <<
                                             this->parameters->pool1_id << "_N" << "\t" << this->parameters->pool1_id << "_I" << "\t" <<
                                             this->parameters->pool2_id << "_A" << "\t" << this->parameters->pool2_id << "_T" << "\t" <<
                                             this->parameters->pool2_id << "_C" << "\t" << this->parameters->pool2_id << "_G" << "\t" <<
                                             this->parameters->pool2_id << "_N" << "\t" << this->parameters->pool2_id << "_I" << "\n";
    }

    if (this->parameters->output_snps_win) {
        this->snps_window_output_file.open(this->parameters->output_prefix);
        this->snps_window_output_file.file << "Contig" << "\t" << "Position" << "\t" << this->parameters->pool1_id << "\t" << this->parameters->pool2_id << "\n";
    }

    if (this->parameters->output_depth) {
        this->depth_output_file.open(this->parameters->output_prefix);
        this->depth_output_file.file << "Contig" << "\t" << "Position" << "\t" << this->parameters->pool1_id << "_abs" << "\t" << this->parameters->pool2_id << "_abs" << "\t"
                                     << this->parameters->pool1_id << "_rel" << "\t" << this->parameters->pool2_id << "_rel" <<"\n";
    }

    if (this->parameters->output_genes) {
        this->genes_output_file.open(this->parameters->output_prefix);
        this->genes_output_file.file << "Contig" << "\t" << "Start" << "\t" << "End" << "\t" << "ID" << "\t" << "Name" << "\t" << "Product" << "\t" <<
                                        this->parameters->pool1_id << "_depth" << "\t" << this->parameters->pool1_id << "_depth_corr" << "\t" <<
                                        this->parameters->pool1_id << "_depth_coding" << "\t" << this->parameters->pool1_id << "_depth_coding_corr" << "\t" <<
                                        this->parameters->pool1_id << "_depth_noncoding" << "\t" << this->parameters->pool1_id << "_depth_noncoding_corr" << "\t" <<
                                        this->parameters->pool2_id << "_depth" << "\t" << this->parameters->pool2_id << "_depth_corr" << "\t" <<
                                        this->parameters->pool2_id << "_depth_coding" << "\t" << this->parameters->pool2_id << "_depth_coding_corr" << "\t" <<
                                        this->parameters->pool2_id << "_depth_noncoding" << "\t" << this->parameters->pool2_id << "_depth_noncoding_corr" << "\t" <<
                                        this->parameters->pool1_id << "_snps" << "\t" << this->parameters->pool1_id << "_snps_coding" << "\t" <<
                                        this->parameters->pool1_id << "_snps_noncoding" << "\t" << this->parameters->pool2_id << "_snps" << "\t" <<
                                        this->parameters->pool2_id << "_snps_coding" << "\t" << this->parameters->pool2_id << "_snps_noncoding" << "\n";
    }
}



// Write SNP and nucleotide information if current base is a sex-specific SNP
void OutputFile::open(const std::string& prefix) {

    this->path = prefix + this->suffix + ".tsv";
    this->file.open(this->path);

    if (not this->file.is_open()) {

        std::cerr << "Error: cannot open output file <" << this->path << ">." << std::endl;
        log("Error: cannot open output file <" + this->path + ">.");
        exit(1);

    }

    log("Created output file <" + this->path + ">.");
}



// Write Fst information if current base has fst higher than specified threshold
void OutputHandler::output_fst_position(float fst) {

    this->fst_position_output_file.file << std::fixed << std::setprecision(4) << this->input_data->contig << "\t" << this->input_data->position << "\t" << fst <<  "\n";
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
void OutputHandler::output_snp_position(std::string& pool_id) {

    this->snps_position_output_file.file << this->input_data->contig << "\t" << this->input_data->position << "\t" << pool_id << "\t"
                                         << this->pair_data->pool1.output_frequencies() << "\t"
                                         << this->pair_data->pool2.output_frequencies() << "\n";
}



// Write SNP and nucleotide information for the current window
void OutputHandler::output_snp_window(uint32_t snps_total[2]) {

    this->snps_window_output_file.file << this->input_data->contig << "\t" << this->input_data->position - this->parameters->window_range << "\t"
                                       << snps_total[0] << "\t"
                                       << snps_total[1] << "\n";
}


// Write depth information at the end of the analysis
void OutputHandler::output_depth(float* average_depth) {

    log("Depth data output started.");

    float window_size = 0;

    for (auto const& contig : *this->depth) {

        for (auto const& position: contig.second) {

            window_size = position.second[2];

            this->depth_output_file.file << contig.first << "\t" << position.first << "\t"
                                         << float(position.second[0] / window_size) << "\t"
                                         << float(position.second[1] / window_size) << "\t"
                                         << std::fixed << std::setprecision(2)
                                         << float((position.second[0] / window_size)/ average_depth[0]) << "\t"
                                         << float((position.second[1] / window_size)/ average_depth[1]) << "\n";

        }
    }

    log("Depth data output ended without errors.");
}


// Write genes information at the end of the analysis
void OutputHandler::output_genes(float* average_depth) {

    log("Genes data output started.");

    float depth_correction_males = (average_depth[0] + average_depth[1]) / 2 / average_depth[0];
    float depth_correction_females = (average_depth[0] + average_depth[1]) / 2 / average_depth[1];

    uint gene_length = 0, male_depth = 0, female_depth = 0;

    for (auto gene: *this->genes) {

        gene_length = uint(std::stoi(gene.second.end) -  std::stoi(gene.second.start));
        male_depth = (gene.second.depth[4 + 0]) / gene_length;
        female_depth = (gene.second.depth[4 + 1]) / gene_length;
        gene.second.noncoding_length = gene_length - gene.second.coding_length;

        if (gene.second.coding_length == 0) {

            gene.second.depth[2 * 0] = 0;
            gene.second.depth[2 * 1] = 0;

        } else {

            (gene.second.noncoding_length > 0) ? gene.second.depth[2 * 0] /= gene.second.noncoding_length : gene.second.depth[2 * 0] = 0;
            (gene.second.noncoding_length > 0) ? gene.second.depth[2 * 1] /= gene.second.noncoding_length: gene.second.depth[2 * 1] = 0;

        }

        if (gene.second.noncoding_length == 0) {

            gene.second.depth[2 * 0 + 1] = 0;
            gene.second.depth[2 * 1 + 1] = 0;

        } else {

            (gene.second.coding_length > 0) ? gene.second.depth[2 * 0 + 1] /= gene.second.coding_length : gene.second.depth[2 * 0 + 1] = 0;
            (gene.second.coding_length > 0) ? gene.second.depth[2 * 1 + 1] /= gene.second.coding_length: gene.second.depth[2 * 1 + 1] = 0;

        }

        this->genes_output_file.file << gene.second.contig << "\t" << gene.second.start << "\t" << gene.second.end << "\t"
                                     <<gene.second.id << "\t" << gene.second.name << "\t" << gene.second.product << "\t"
                                     << male_depth << "\t" << int(male_depth * depth_correction_males) << "\t"
                                     << gene.second.depth[2 * 0 + 1] << "\t" << int(gene.second.depth[2 * 0 + 1] * depth_correction_males) << "\t"
                                     << gene.second.depth[2 * 0] << "\t" << int(gene.second.depth[2 * 0] * depth_correction_males) << "\t"
                                     << female_depth << "\t" << int(female_depth * depth_correction_females) << "\t"
                                     << gene.second.depth[2 * 1 + 1] << "\t" << int(gene.second.depth[2 * 1 + 1] * depth_correction_females) << "\t"
                                     << gene.second.depth[2 * 1] << "\t" << int(gene.second.depth[2 * 1] * depth_correction_females) << "\t"
                                     << gene.second.snps[4 + 0] << "\t" << gene.second.snps[2 * 0 + 1] << "\t" << gene.second.snps[2 * 0] << "\t"
                                     << gene.second.snps[4 + 1] << "\t" << gene.second.snps[2 * 1 + 1] << "\t" << gene.second.snps[2 * 1] << "\n";
    }

    log("Genes data output ended without errors.");

}
