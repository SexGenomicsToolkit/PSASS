#include "arg_parser.h"


ArgParser::ArgParser(int &argc, char **argv) {

    for (auto i=1; i < argc; ++i) this->fields.push_back(std::string(argv[i]));

    if (this->contains("-h")) {
        this->usage();
        exit(0);
    }

    if (!this->contains("-i")){
        std::cout << std::endl << "** Error: no input file specified" << std::endl;
        this->usage();
        exit(0);
    }

    if (!this->contains("-o")){
        std::cout << std::endl << "** Error: no output file specified" << std::endl;
        this->usage();
        exit(0);
    }
}



void ArgParser::set_parameters(Parameters& parameters) {

    parameters.min_depth = std::stoul(this->set_value("-d"));
    if (parameters.min_depth > 0) parameters.min_depth -= 1; //Set value - 1 to avoid doing ">=" later
    parameters.min_fst = std::stof(this->set_value("-f"));
    if (parameters.min_fst > 0) parameters.min_fst -= 0.000001; //Set value - 0.000001 to avoid doing ">=" later
    parameters.snp_range = std::stof(this->set_value("-s"));
    if (parameters.snp_range < 1) parameters.snp_range += 0.000001; //Set value + 0.000001 to avoid doing "<=" later
    parameters.window_size = std::stoul(this->set_value("-w"));
    parameters.output_resolution = std::stoul(this->set_value("-r"));
    parameters.input_file_path = this->set_value("-i");
    parameters.output_file_path = this->set_value("-o");
    parameters.male_pool = std::stoul(this->set_value("-m"));

    parameters.input_file.open(parameters.input_file_path);

    if (not parameters.input_file.is_open()) {
        std::cout << "Error: cannot open input file (" << parameters.input_file_path << ")." << std::endl;
        exit(0);
    }

    if (parameters.output_file_path != "") parameters.output_file_path += "_";
    std::string snps_output_file_path = parameters.output_file_path + "snps.tsv";
    parameters.snps_output_file.open(snps_output_file_path);
    if (not parameters.snps_output_file.is_open()) {
        std::cout << "Error: cannot open SNPs output file (" << snps_output_file_path << ")." << std::endl;
        exit(0);
    }

    std::string fst_threshold_output_file_path = parameters.output_file_path + "fst_threshold.tsv";
    parameters.fst_threshold_output_file.open(fst_threshold_output_file_path);
    if (not parameters.fst_threshold_output_file.is_open()) {
        std::cout << "Error: cannot open Fst threshold output file (" << fst_threshold_output_file_path << ")." << std::endl;
        exit(0);
    }

    std::string fst_window_output_file_path = parameters.output_file_path + "fst_window.tsv";
    parameters.fst_window_output_file.open(fst_window_output_file_path);
    if (not parameters.fst_window_output_file.is_open()) {
        std::cout << "Error: cannot open Fst window output file (" << fst_window_output_file_path << ")." << std::endl;
        exit(0);
    }

    parameters.fst_threshold_output_file << "Contig" << "\t" << "Position" << "\t" << "Fst" << "\n";
    parameters.fst_window_output_file << "Contig" << "\t" << "Position" << "\t" << "Fst" << "\n";
    parameters.snps_output_file << "Contig" << "\t" << "Position" << "\t" << "Males" << "\t" << "Females" << "\n";

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

    std::cout << std::endl << "Usage: poolsex [options] -i input_file.sync -o output_file" << std::endl;
    std::cout << std::endl << "Options:" << std::endl << std::endl;
    for (auto o: this->options) std::cout << "\t" << o.first << " <" << o.second[1] << ">  " << o.second[2] << "  [" << o.second[0] << "]" << std::endl;
    std::cout << std::endl;
}


void ArgParser::print_parameters() {

    std::cout << "\n- Parameters:\n";
    for (auto o: this->options) {
        if (o.first != "-h") std::cout << "\t" << "- " << o.second[2] << " : " << this->set_value(o.first) << "\n";
    }

    std::cout << "\n";
}
