#pragma once
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "parameters.h"
#include "utils.h"


class ArgParser {

    public:

        // Options: flag -> [default, type, help message]
        std::map<std::string, std::vector<std::string>> options { {"--help", {"0", "bool", "Prints this message"} },

                                                                  {"--output-fst-pos", {"1", "bool", "Output fst positions"} },
                                                                  {"--output-fst-win", {"1", "bool", "Output fst sliding window"} },
                                                                  {"--output-snps-pos", {"1", "bool", "Output snps positions"} },
                                                                  {"--output-snps-win", {"1", "bool", "Output snps sliding window"} },
                                                                  {"--output-depth", {"1", "bool", "Output depth for each pool"} },
                                                                  {"--male-pool", {"2", "int", "Male pool (1/2)"} } ,

                                                                  {"--min-depth", {"10", "int", "Minimum depth to consider a site"} },
                                                                  {"--min-fst", {"0.25", "float", "FST threshold"} },
                                                                  {"--freq-het", {"0.5", "float", "Frequency of a sex-linked SNP in the heterogametic sex"} },
                                                                  {"--freq-hom", {"1", "float", "Frequency of a sex-linked SNP in the homogametic sex"} },
                                                                  {"--range-het", {"0.1", "float", "Range of frequency for a sex-linked SNP in the heterogametic sex"} },
                                                                  {"--range-hom", {"0.05", "float", "Range of frequency for a sex-linked SNP in the homogametic sex"} },
                                                                  {"--window-size", {"100000", "int", "Size of the sliding window (in bp)"} },
                                                                  {"--output-resolution", {"10000", "int", "Output resolution (in bp)"} } ,
                                                                  {"--group-snps", {"1", "bool", "Group consecutive snps to count them as a single polymorphism"} },

                                                                  {"--input-file", {"", "string", "Input file (psass convert output or popoolation sync file)"} },
                                                                  {"--input-format", {"", "string", "Input format (psass/popoolation)"} },
                                                                  {"--output-prefix", {"", "string", "Full prefix (including path) for output files"} },
                                                                  {"--gff-file", {"", "string", "GFF file for gene-specific output"} }
                                                                };

        std::vector<std::string> print_order {"#Input", "--input-file", "--input-format", "--gff-file", "--male-pool", "#Output", "--output-prefix", "--output-fst-pos",
                                              "--output-fst-win", "--output-snps-pos", "--output-snps-win", "--output-depth", "#Computations", "--min-depth", "--min-fst",
                                              "--freq-het", "--range-het", "--freq-hom", "--range-hom", "--window-size", "--output-resolution"};

        ArgParser(int &argc, char **argv);
        void set_parameters(Parameters& parameters);
        const std::string get_value(const std::string& setting) const;
        bool contains(const std::string &option) const ;
        const std::string set_value(const std::string& field);
        void usage();
        void print_parameters();
        void output_parameters(std::ofstream& output_file);

    private:

        std::vector<std::string> fields;
};
