#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include "parameters.h"


class ArgParser {

    public:

        // Options: flag -> [default, type, help message]
        std::map<std::string, std::vector<std::string>> const options { {"-h", {"0", "bool", "Prints this message"} },
                                                                        {"-d", {"10", "int", "Minimum depth to consider a site"} },
                                                                        {"-f", {"0.2", "float", "FST threshold"} },
                                                                        {"-s", {"0.1", "float", "Range to consider a sex-linked SNP"} },
                                                                        {"-x", {"0.1", "float", "Range to consider a fixed SNP"} },
                                                                        {"-w", {"100000", "int", "Size of the sliding window (in bp)"} },
                                                                        {"-r", {"500", "int", "Output resolution (in bp)"} } ,
                                                                        {"-m", {"1", "int", "Male pool (1/2)"} } ,
                                                                        {"-i", {"", "string", "Input file (popoolation sync file)"} },
                                                                        {"-o", {"", "string", "Output file"} },
                                                                        {"-c", {"1", "bool", "Output coverage"} },
                                                                        {"-p", {"1", "bool", "Output snps positions"} }
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
