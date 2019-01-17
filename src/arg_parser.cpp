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



void ArgParser::set_parameters(Parameters& parameters) {

    parameters.input_file_path = this->set_value("--input-file");
    parameters.output_prefix = this->set_value("--output-prefix");


    // Retrieve parameter values
    parameters.min_depth = uint(std::stoul(this->set_value("--min-depth")));
    parameters.min_fst = std::stof(this->set_value("--min-fst"));
    parameters.freq_het = std::stof(this->set_value("--freq-het"));
    parameters.freq_hom = std::stof(this->set_value("--freq-hom"));
    parameters.range_het = std::stof(this->set_value("--range-het"));
    parameters.range_hom = std::stof(this->set_value("--range-hom"));
    parameters.min_het = parameters.freq_het - parameters.range_het;
    parameters.max_het = parameters.freq_het + parameters.range_het;
    parameters.min_hom = parameters.freq_hom- parameters.range_hom;
    parameters.group_snps = this->set_value("--group-snps") != "0";
    parameters.output_resolution = uint(std::stoul(this->set_value("--output-resolution")));
    parameters.window_size = uint(std::stoul(this->set_value("--window-size")));
    parameters.window_range = parameters.window_size / 2;
    parameters.male_pool = uint(std::stoul(this->set_value("--male-pool")));
    parameters.output_fst_pos = this->set_value("--output-fst-pos") != "0";
    parameters.output_fst_win = this->set_value("--output-fst-win") != "0";
    parameters.output_snps_pos = this->set_value("--output-snps-pos") != "0";
    parameters.output_snps_win = this->set_value("--output-snps-win") != "0";
    parameters.output_depth = this->set_value("--output-depth") != "0";

    // Open input file
    parameters.input_file.open(parameters.input_file_path);
    if (not parameters.input_file.is_open()) {
        std::cerr << "Error: cannot open input file (" << parameters.input_file_path << ")." << std::endl;
        exit(1);
    }


    // Open GFF file
    if (this->contains("--gff-file")) {
        parameters.output_genes = true;
        parameters.gff_file_path = this->set_value("--gff-file");
        parameters.gff_file.open(parameters.gff_file_path);
        if (not parameters.gff_file.is_open()) {
            std::cerr << "Error: cannot open gff file (" << parameters.gff_file_path << ")." << std::endl;
            exit(1);
        }
    }
}



void ArgParser::usage() {

    std::cout << std::endl << "Usage: poolsex [options] --input-file input_file.sync --output-prefix output_prefix" << std::endl;
    std::cout << std::endl << "Options:" << std::endl << std::endl;

    std::cout << "## Input / output " << std::endl;
    std::cout << "--input-file           <string>    Input file (popoolation sync file)                                    [\"\"]" << std::endl;
    std::cout << "--output-prefix        <string>    Full prefix (including path) for output files                         [\"\"]" << std::endl;
    std::cout << "--gff-file             <string>    GFF file for gene-specific output                                     [\"\"]" << std::endl;
    std::cout << "--output-fst-pos       <bool>      If true, output fst positions (0/1)                                   [1]" << std::endl;
    std::cout << "--output-fst-win       <bool>      If true, output fst sliding window (0/1)                              [1]" << std::endl;
    std::cout << "--output-snps-pos      <bool>      If true, output snps positions (0/1)                                  [1]" << std::endl;
    std::cout << "--output-snps-win      <bool>      If true, output snps sliding window (0/1)                             [1]" << std::endl;
    std::cout << "--output-depth         <bool>      If true, output depth(0/1)                                            [1]" << std::endl;
    std::cout << "--male-pool            <int>       Male pool (1/2)                                                       [2]" << std::endl << std::endl;
    std::cout << "## Analysis" << std::endl;
    std::cout << "--min-depth            <int>       Minimum depth to consider a site                                      [10]" << std::endl;
    std::cout << "--min-fst              <float>     FST threshold                                                         [0.25]" << std::endl;
    std::cout << "--freq-het             <float>     Frequency of a sex-linked SNP in the heterogametic sex                [0.5]" << std::endl;
    std::cout << "--freq-hom             <float>     Frequency of a sex-linked SNP in the homogametic sex                  [1]" << std::endl;
    std::cout << "--range-het            <float>     Range of frequency for a sex-linked SNP in the heterogametic sex      [0.1]" << std::endl;
    std::cout << "--range-hom            <float>     Range of frequency for a sex-linked SNP in the homogametic sex        [0.05]" << std::endl;
    std::cout << "--window-size          <int>       Size of the sliding window (in bp)                                    [100000]" << std::endl;
    std::cout << "--output-resolution    <int>       Output resolution (in bp)                                             [500]" << std::endl;
    std::cout << "--group-snps           <bool>      Group consecutive snps to count them as a single polymorphism (0/1)   [0]" << std::endl;
    std::cout << std::endl;
}



void ArgParser::print_parameters() {

    std::cout << "\n- Parameters:\n";
    for (auto o: this->options) {
        if (o.first != "--help") std::cout << "\t" << "- " << o.second[2] << " : " << this->set_value(o.first) << "\n";
    }

    std::cout << "\n";
}



void ArgParser::output_parameters(std::ofstream& output_file) {

    output_file << "\nPSASS parameters:\n";

    for (std::string o: this->print_order) {
        if (o.substr(0,1) == "#") output_file << o << "\n";
        else if (o != "--help") output_file << " - " << this->options[o][2] << " : " << this->set_value(o) << "\n";
    }

    output_file << std::endl;
}
