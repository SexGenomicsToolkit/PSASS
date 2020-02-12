#include "output_handler.h"

OutputHandler::OutputHandler(Parameters& parameters) {

    // Open output file objects
    if (parameters.fst_pos_file_path != "") {
        this->open_output_file(this->fst_position_output_file, parameters.fst_pos_file_path);
        this->fst_position_output_file << "Contig" << "\t" << "Position" << "\t" << "Length" << "\t" << "Fst" << "\n";
    }

    if (parameters.snp_pos_file_path != "") {
        this->open_output_file(this->snp_position_output_file, parameters.snp_pos_file_path);
        this->snp_position_output_file << "Contig" << "\t" << "Position" << "\t" << "Length" << "\t" << "Pool" << "\t" <<
                                          parameters.pool1_id << "_A" << "\t" << parameters.pool1_id << "_T" << "\t" <<
                                          parameters.pool1_id << "_C" << "\t" << parameters.pool1_id << "_G" << "\t" <<
                                          parameters.pool1_id << "_N" << "\t" << parameters.pool1_id << "_O" << "\t" <<
                                          parameters.pool2_id << "_A" << "\t" << parameters.pool2_id << "_T" << "\t" <<
                                          parameters.pool2_id << "_C" << "\t" << parameters.pool2_id << "_G" << "\t" <<
                                          parameters.pool2_id << "_N" << "\t" << parameters.pool2_id << "_O" << "\n";
    }

    if (parameters.genes_file_path != "") {
        this->open_output_file(this->genes_output_file, parameters.genes_file_path);
        this->genes_output_file << "Contig" << "\t" << "Start" << "\t" << "End" << "\t" << "ID" << "\t" << "Name" << "\t" << "Product" << "\t" <<
                                   parameters.pool1_id << "_depth" << "\t" << parameters.pool1_id << "_depth_corr" << "\t" <<
                                   parameters.pool1_id << "_depth_coding" << "\t" << parameters.pool1_id << "_depth_coding_corr" << "\t" <<
                                   parameters.pool1_id << "_depth_noncoding" << "\t" << parameters.pool1_id << "_depth_noncoding_corr" << "\t" <<
                                   parameters.pool2_id << "_depth" << "\t" << parameters.pool2_id << "_depth_corr" << "\t" <<
                                   parameters.pool2_id << "_depth_coding" << "\t" << parameters.pool2_id << "_depth_coding_corr" << "\t" <<
                                   parameters.pool2_id << "_depth_noncoding" << "\t" << parameters.pool2_id << "_depth_noncoding_corr" << "\t" <<
                                   parameters.pool1_id << "_snps" << "\t" << parameters.pool1_id << "_snps_coding" << "\t" <<
                                   parameters.pool1_id << "_snps_noncoding" << "\t" << parameters.pool2_id << "_snps" << "\t" <<
                                   parameters.pool2_id << "_snps_coding" << "\t" << parameters.pool2_id << "_snps_noncoding" << "\n";
    }

    this->open_output_file(this->window_output_file, parameters.output_file_path);
    this->window_output_file << "Contig" << "\t" << "Position" << "\t" << "Length" << "\t"
                             << "Snps_" << parameters.pool1_id << "\t" << "Snps_" << parameters.pool2_id << "\t"
                             << "Fst" << "\t"
                             << "Abs_depth_" << parameters.pool1_id << "\t" << "Abs_depth_" << parameters.pool2_id << "\t"
                             << "Rel_depth_" << parameters.pool1_id << "\t" << "Rel_depth_" << parameters.pool2_id << "\t"
                             << "Depth_ratio" << "\n";

    this->pool_id[0] = parameters.pool1_id;
    this->pool_id[1] = parameters.pool2_id;
    this->min_depth = parameters.min_depth;
}


void OutputHandler::open_output_file(std::ofstream& output_file, std::string& path) {

    output_file.open(path);

    if (not output_file.is_open()) {

        std::cerr << "Error: cannot open output file <" << path << ">." << std::endl;
        log("Error: cannot open output file <" + path + ">.");
        exit(1);

    }

    log("Created output file <" + path + ">.");
}



// Output single-base information for snp and fst output files
void OutputHandler::output_bases(std::vector<PointOutputData> base_output_data, std::unordered_map<std::string, uint64_t>& contig_lengths) {

    for (auto base: base_output_data) {

        if (base.type == 0) {

            this->fst_position_output_file << base.contig << "\t" << base.position << "\t" << contig_lengths[base.contig] << "\t" << std::fixed << std::setprecision(4) << base.base_data.fst << "\n";

        } else {

            this->snp_position_output_file << base.contig << "\t" << base.position << "\t" << contig_lengths[base.contig] << "\t" << this->pool_id[base.type - 1] << "\t"
                                           << base.base_data.pool1.output_frequencies() << "\t"
                                           << base.base_data.pool2.output_frequencies() << "\n";

        }
    }
}



// Write genes information at the end of the analysis
void OutputHandler::output_genes(std::unordered_map<std::string, Gene>& genes, float* average_depth) {

    log("Genes data output started.");

    float depth_correction_males = (average_depth[0] + average_depth[1]) / 2 / average_depth[0];
    float depth_correction_females = (average_depth[0] + average_depth[1]) / 2 / average_depth[1];

    uint gene_length = 0, male_depth = 0, female_depth = 0;

    for (auto gene: genes) {

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

        this->genes_output_file << gene.second.contig << "\t" << gene.second.start << "\t" << gene.second.end << "\t"
                                << gene.second.id << "\t" << gene.second.name << "\t" << gene.second.product << "\t"
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


void OutputHandler::output_window(std::map<std::string, std::map<uint, float[6]>>& output_data, float* average_depth, std::unordered_map<std::string, uint64_t>& contig_lengths) {

    for (auto const& contig : output_data) {
        for (auto const& position: contig.second) {
            this->window_output_file << contig.first << "\t" << position.first << "\t" << contig_lengths[contig.first] << "\t"
                                     << std::fixed << std::setprecision(0)
                                     << uint(position.second[3]) << "\t" << uint(position.second[4]) << "\t"
                                     << std::fixed << std::setprecision(4)
                                     << position.second[5] << "\t"
                                     << std::fixed << std::setprecision(0)
                                     << float(position.second[0] / position.second[2]) << "\t"
                                     << float(position.second[1] / position.second[2]) << "\t"
                                     << std::fixed << std::setprecision(2)
                                     << float((position.second[0] / position.second[2])/ average_depth[0]) << "\t"
                                     << float((position.second[1] / position.second[2])/ average_depth[1]) << "\t";

            if (position.second[0] >= this->min_depth and position.second[1] >= this->min_depth) {

                this->window_output_file << float(position.second[0] / position.second[1]) << "\n";

            } else {

                this->window_output_file << 1.00 << "\n";

            }
        }
    }
}
