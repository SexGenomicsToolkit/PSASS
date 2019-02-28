#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include "logs.h"
#include "utils.h"

class PileupConverter {

    public:

        PileupConverter(int argc, char *argv[]);
        void run();

    private:

        std::unordered_map<char, uint> bases {{'A', 0}, {'T', 1}, {'G', 2}, {'C', 3}, {'N', 4}, {'I', 5},
                                              {'a', 0}, {'t', 1}, {'g', 2}, {'c', 3}, {'n', 4}};

        // Efficient file reading parameters
        char buff[2048];
        int buff_size = sizeof(this->buff);
        long k = 0;

        // I/O parameters
        bool from_stdin = false;
        std::ifstream input_file;
        std::ofstream output_file;

        // Syncfile parsing : current field and subfield in a line, and current (sub)field value is stored in temp
        uint field = 0;
        std::string temp = "";
        uint n_lines = 0;

        // Fields
        std::string contig = "";
        uint position = 0;
        char ref_allele = 'N';
        uint depth[2] = {0, 0};
        uint pool[2][6] {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};

        // Base processing
        bool read_begin = false;
        bool next_indel = false;
        std::string tmp_indel_size = "";
        uint remaining_indel = 0;

        void process_line();
        void process_field();
        void process_base(char base, bool pool);
        void usage();
};
