#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include "parameters.h"
#include "utils.h"


class ArgParser {

    public:

        // Options: flag -> [default, type, help message]
        std::map<std::string, std::vector<std::string>> const options { {"--help", {"0", "bool", "Prints this message"} },

                                                                        {"--output-fst-pos", {"1", "bool", "Output fst positions"} },
                                                                        {"--output-fst-win", {"1", "bool", "Output fst sliding window"} },
                                                                        {"--output-snps-pos", {"1", "bool", "Output snps positions"} },
                                                                        {"--output-snps-win", {"1", "bool", "Output snps sliding window"} },
                                                                        {"--output-coverage", {"1", "bool", "Output coverage"} },
                                                                        {"--male-pool", {"2", "int", "Male pool (1/2)"} } ,

                                                                        {"--min-depth", {"10", "int", "Minimum depth to consider a site"} },
                                                                        {"--min-fst", {"0.25", "float", "FST threshold"} },
                                                                        {"--freq-het", {"0.5", "float", "Frequency of a sex-linked SNP in the heterogametic sex"} },
                                                                        {"--freq-hom", {"1", "float", "Frequency of a sex-linked SNP in the homogametic sex"} },
                                                                        {"--range-het", {"0.1", "float", "Range of frequency for a sex-linked SNP in the heterogametic sex"} },
                                                                        {"--range-hom", {"0.05", "float", "Range of frequency for a sex-linked SNP in the homogametic sex"} },
                                                                        {"--window-size", {"100000", "int", "Size of the sliding window (in bp)"} },
                                                                        {"--output-resolution", {"500", "int", "Output resolution (in bp)"} } ,

                                                                        {"--input-file", {"", "string", "Input file (popoolation sync file)"} },
                                                                        {"--output-prefix", {"", "string", "Full prefix (including path) for output files"} },
                                                                        {"--gff-file", {"", "string", "GFF file for gene-specific output"} }
                                                                      };

        ArgParser(int &argc, char **argv);

        void set_parameters(Parameters& parameters);

        const std::string get_value(const std::string& setting) const;

        bool contains(const std::string &option) const ;

        const std::string set_value(const std::string& field);

        void usage();

        void print_parameters();

    private:

        std::vector<std::string> fields;
};
