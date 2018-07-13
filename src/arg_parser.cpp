#include "arg_parser.h"


ArgParser::ArgParser(int &argc, char **argv) {

    for (auto i=1; i < argc; ++i) this->fields.push_back(std::string(argv[i]));

    if (this->contains("--help")) {
        this->usage();
        exit(1);
    }

    if (!this->contains("--input-file")) {
        std::cout << std::endl << "** Error: no input file specified" << std::endl;
        this->usage();
        exit(1);
    }

    if (!this->contains("--output-prefix")) {
        std::cout << std::endl << "** Error: no output prefix specified" << std::endl;
        this->usage();
        exit(1);
    }
}



void ArgParser::set_parameters(Parameters& parameters) {

    parameters.input_file_path = this->set_value("--input-file");
    parameters.output_prefix = this->set_value("--output-prefix");
    parameters.log_file.open(parameters.output_prefix + ".log");

    write_log("Retrieving parameter values ... \n", parameters.log_file, true, true);

    // Retrieve parameter values
    parameters.min_depth = std::stoul(this->set_value("--min-depth"));
    write_log(" - Min depth: ", parameters.log_file, false, false);
    write_log(parameters.min_depth, parameters.log_file, false, true);

    parameters.min_fst = std::stof(this->set_value("--min-fst"));
    write_log(" - Min FST: ", parameters.log_file, false, false);
    write_log(parameters.min_fst, parameters.log_file, false, true);

    parameters.freq_het = std::stof(this->set_value("--freq-het"));
    write_log(" - Heterozygous frequency: ", parameters.log_file, false, false);
    write_log(parameters.freq_het, parameters.log_file, false, true);

    parameters.freq_hom = std::stof(this->set_value("--freq-hom"));
    write_log(" - Homozygous frequency: ", parameters.log_file, false, false);
    write_log(parameters.freq_hom, parameters.log_file, false, true);

    parameters.range_het = std::stof(this->set_value("--range-het"));
    write_log(" - Heterozygous range: ", parameters.log_file, false, false);
    write_log(parameters.range_het, parameters.log_file, false, true);

    parameters.range_hom = std::stof(this->set_value("--range-hom"));
    write_log(" - Homozygous range: ", parameters.log_file, false, false);
    write_log(parameters.range_hom, parameters.log_file, false, true);

    parameters.window_size = std::stoul(this->set_value("--window-size"));
    write_log(" - Window size: ", parameters.log_file, false, false);
    write_log(parameters.window_size, parameters.log_file, false, true);

    parameters.output_resolution = std::stoul(this->set_value("--output-resolution"));
    write_log(" - Output resolution: ", parameters.log_file, false, false);
    write_log(parameters.output_resolution, parameters.log_file, false, true);

    parameters.male_pool = std::stoul(this->set_value("--male-pool"));
    write_log(" - Male pool: ", parameters.log_file, false, false);
    write_log(parameters.male_pool, parameters.log_file, false, true);

    parameters.output_fst_pos = this->set_value("--output-fst-pos") != "0";
    write_log(" - Output Fst positions: ", parameters.log_file, false, false);
    write_log(parameters.output_fst_pos, parameters.log_file, false, true);

    parameters.output_fst_win = this->set_value("--output-fst-win") != "0";
    write_log(" - Output Fst window: ", parameters.log_file, false, false);
    write_log(parameters.output_fst_win, parameters.log_file, false, true);

    parameters.output_snps_pos = this->set_value("--output-snps-pos") != "0";
    write_log(" - Output Snps positions: ", parameters.log_file, false, false);
    write_log(parameters.output_snps_pos, parameters.log_file, false, true);

    parameters.output_snps_win = this->set_value("--output-snps-win") != "0";
    write_log(" - Output Snps window: ", parameters.log_file, false, false);
    write_log(parameters.output_snps_win, parameters.log_file, false, true);

    parameters.output_coverage = this->set_value("--output-coverage") != "0";
    write_log(" - Output coverage: ", parameters.log_file, false, false);
    write_log(parameters.output_coverage, parameters.log_file, false, true);

    // Open input file
    write_log("\n", parameters.log_file, false, false);
    write_log("Opening input sync file: ", parameters.log_file, true, false);
    parameters.input_file.open(parameters.input_file_path);
    if (not parameters.input_file.is_open()) {
        std::cout << "Error: cannot open input file (" << parameters.input_file_path << ")." << std::endl;
        write_log("cannot open input file (" + parameters.input_file_path + ").", parameters.log_file, false, true);
        exit(1);
    }
    write_log("OK", parameters.log_file, false, true);


    // Open GFF file
    if (this->contains("--gff-file")) {
        write_log("Opening gff file: ", parameters.log_file, true, false);
        parameters.gff_file_path = this->set_value("--gff-file");
        parameters.gff_file.open(parameters.gff_file_path);
        if (not parameters.gff_file.is_open()) {
            std::cout << "Error: cannot open gff file (" << parameters.gff_file_path << ")." << std::endl;
            write_log("cannot open input file (" + parameters.gff_file_path + ").", parameters.log_file, false, true);
            exit(1);
        }
        write_log("OK", parameters.log_file, false, true);
    }


    // Create base output file path
    if (parameters.output_prefix != "") parameters.output_prefix += "_";

    write_log("Creating basic output files ... \n", parameters.log_file, true, true);

    // Position Fst output file
    if (parameters.output_fst_pos) {
        write_log(" - FST position file: ", parameters.log_file, false, false);
        std::string fst_position_output_file_path = parameters.output_prefix + "position_fst.tsv";
        parameters.fst_pos_output_file.open(fst_position_output_file_path);
        if (not parameters.fst_pos_output_file.is_open()) {
            std::cout << "Error: cannot open Fst threshold output file (" << fst_position_output_file_path << ")." << std::endl;
            write_log("cannot open input file (" + fst_position_output_file_path + ").", parameters.log_file, false, true);
            exit(1);
        }
        parameters.fst_pos_output_file << "Contig" << "\t" << "Position" << "\t" << "Fst" << "\n";
        write_log("OK", parameters.log_file, false, true);
    }

    // Window Fst output file
    if (parameters.output_fst_win) {
        write_log(" - FST window file: ", parameters.log_file, false, false);
        std::string fst_win_output_file_path = parameters.output_prefix + "window_fst.tsv";
        parameters.fst_win_output_file.open(fst_win_output_file_path);
        if (not parameters.fst_win_output_file.is_open()) {
            std::cout << "Error: cannot open Fst window output file (" << fst_win_output_file_path << ")." << std::endl;
            write_log("cannot open input file (" + fst_win_output_file_path + ").", parameters.log_file, false, true);
            exit(1);
        }
        parameters.fst_win_output_file << "Contig" << "\t" << "Position" << "\t" << "Fst" << "\n";
        write_log("OK", parameters.log_file, false, true);
    }

    // Position SNPs output file
    if (parameters.output_snps_pos) {
        write_log(" - SNPs position file: ", parameters.log_file, false, false);
        std::string snps_pos_output_file_path = parameters.output_prefix + "position_snp.tsv";
        parameters.snps_pos_output_file.open(snps_pos_output_file_path);
        if (not parameters.snps_pos_output_file.is_open()) {
            std::cout << "Error: cannot open snps pos output file (" << snps_pos_output_file_path << ")." << std::endl;
            write_log("cannot open input file (" + snps_pos_output_file_path + ").", parameters.log_file, false, true);
            exit(1);
        }
        parameters.snps_pos_output_file << "Contig" << "\t" << "Position" << "\t" << "Sex" << "\t" <<
                                           "M_A" << "\t" << "M_T" << "\t" << "M_G" << "\t" << "M_C" << "\t" << "M_I" << "\t" <<
                                           "F_A" << "\t" << "F_T" << "\t" << "F_G" << "\t" << "F_C" << "\t" << "F_I" << "\n";
        write_log("OK", parameters.log_file, false, true);
    }

    // Window SNPs output file
    if (parameters.output_snps_win) {
        write_log(" - SNPs window file: ", parameters.log_file, false, false);
        std::string snps_win_output_file_path = parameters.output_prefix + "window_snp.tsv";
        parameters.snps_win_output_file.open(snps_win_output_file_path);
        if (not parameters.snps_win_output_file.is_open()) {
            std::cout << "Error: cannot open SNPs output file (" << snps_win_output_file_path << ")." << std::endl;
            write_log("cannot open input file (" + snps_win_output_file_path + ").", parameters.log_file, false, true);
            exit(1);
        }
        parameters.snps_win_output_file << "Contig" << "\t" << "Position" << "\t" << "Males" << "\t" << "Females" << "\n";
        write_log("OK", parameters.log_file, false, true);
    }

    // Coverage output file
    if (parameters.output_coverage) {
        write_log(" - Coverage file: ", parameters.log_file, false, false);
        std::string coverage_output_file_path = parameters.output_prefix + "coverage.tsv";
        parameters.coverage_output_file.open(coverage_output_file_path);
        if (not parameters.coverage_output_file.is_open()) {
            std::cout << "Error: cannot open coverage output file (" << coverage_output_file_path << ")." << std::endl;
            write_log("cannot open input file (" + coverage_output_file_path + ").", parameters.log_file, false, true);
            exit(1);
        }
        parameters.coverage_output_file << "Contig" << "\t" << "Position" << "\t" << "Males_rel" << "\t" << "Females_rel" << "\t" << "Males_abs" << "\t" << "Females_abs" << "\n";
        write_log("OK", parameters.log_file, false, true);
    }


    // Genes output file
    if (this->contains("--gff-file")) {
        write_log(" - Genes file: ", parameters.log_file, false, false);
        std::string genes_output_file_path = parameters.output_prefix + "genes.tsv";
        parameters.genes_output_file.open(genes_output_file_path);
        if (not parameters.genes_output_file.is_open()) {
            std::cout << "Error: cannot open genes output file (" << genes_output_file_path << ")." << std::endl;
            write_log("cannot open input file (" + genes_output_file_path + ").", parameters.log_file, false, true);
            exit(1);
        }
        parameters.genes_output_file << "Contig" << "\t" << "Start" << "\t" << "End" << "\t" <<
                                        "Name" << "\t" << "Product" << "\t" <<
                                        "Cov_males" << "\t" << "Cov_males_coding" << "\t" << "Cov_males_noncoding" << "\t" <<
                                        "Cov_females" << "\t" << "Cov_females_coding" << "\t" << "Cov_females_noncoding" << "\t" <<
                                        "Snp_males" << "\t" << "Snp_males_coding" << "\t" << "Snp_males_noncoding" << "\t" <<
                                        "Snp_females" << "\t" << "Snp_females_coding" << "\t" << "Snp_females_noncoding" << "\n";
        parameters.output_genes = true;
        write_log("OK", parameters.log_file, false, true);
    }
}



const std::string ArgParser::get_value(const std::string& setting) const {

    std::vector<std::string>::const_iterator itr = std::find(this->fields.begin(), this->fields.end(), setting);

    if (itr != this->fields.end() && ++itr != this->fields.end()) {

        return *itr;
    }

    return "";
}



bool ArgParser::contains(const std::string &option) const {

    return std::find(this->fields.begin(), this->fields.end(), option) != this->fields.end();
}



const std::string ArgParser::set_value(const std::string& field) {

    if (this->contains(field)) return this->get_value(field);
    else  return this->options.at(std::string(field))[0];
}



void ArgParser::usage() {

    std::cout << std::endl << "Usage: poolsex [options] --input-file input_file.sync --output-prefix output_prefix" << std::endl;
    std::cout << std::endl << "Options:" << std::endl << std::endl;

    std::cout << "## Input / output " << std::endl;
    std::cout << "--input-file           <string>    Input file (popoolation sync file)                                    [\"\"]" << std::endl;
    std::cout << "--output-prefix        <string>    Full prefix (including path) for output files                         [\"\"]" << std::endl;
    std::cout << "--gff-file             <string>    GFF file for gene-specific output                                     [\"\"]" << std::endl;
    std::cout << "--output-fst-pos       <bool>      Output fst positions                                                  [1]" << std::endl;
    std::cout << "--output-fst-win       <bool>      Output fst sliding window                                             [1]" << std::endl;
    std::cout << "--output-snps-pos      <bool>      Output snps positions                                                 [1]" << std::endl;
    std::cout << "--output-snps-win      <bool>      Output snps sliding window                                            [1]" << std::endl;
    std::cout << "--output-coverage      <bool>      Output coverage                                                       [1]" << std::endl;
    std::cout << "--male-pool            <int>       Male pool (1/2)                                                       [2]" << std::endl << std::endl;;
    std::cout << "## Analysis" << std::endl;
    std::cout << "--min-depth            <int>       Minimum depth to consider a site                                      [10]" << std::endl;
    std::cout << "--min-fst              <float>     FST threshold                                                         [0.25]" << std::endl;
    std::cout << "--freq-het             <float>     Frequency of a sex-linked SNP in the heterogametic sex                [0.5]" << std::endl;
    std::cout << "--freq-hom             <float>     Frequency of a sex-linked SNP in the homogametic sex                  [1]" << std::endl;
    std::cout << "--range-het            <float>     Range of frequency for a sex-linked SNP in the heterogametic sex      [0.1]" << std::endl;
    std::cout << "--range-hom            <float>     Range of frequency for a sex-linked SNP in the homogametic sex        [0.05]" << std::endl;
    std::cout << "--window-size          <int>       Size of the sliding window (in bp)                                    [100000]" << std::endl;
    std::cout << "--output-resolution    <int>       Output resolution (in bp)                                             [500]" << std::endl;

//    for (auto o: this->options) std::cout << "\t" << o.first << " <" << o.second[1] << ">  " << o.second[2] << "  [" << o.second[0] << "]" << std::endl;
    std::cout << std::endl;
}



void ArgParser::print_parameters() {

    std::cout << "\n- Parameters:\n";
    for (auto o: this->options) {
        if (o.first != "--help") std::cout << "\t" << "- " << o.second[2] << " : " << this->set_value(o.first) << "\n";
    }

    std::cout << "\n";
}
