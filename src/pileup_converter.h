#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include "parameters.h"
#include "utils.h"

class PileupConverter {

    public:

        PileupConverter(Parameters& parameters);
        void run();

    private:

        // Efficient file reading parameters
        char buff[8192];
        int buff_size = sizeof(this->buff);
        long k = 0;

        // I/O parameters
        bool from_stdin = false;
        bool to_stdout = true;
        std::ifstream input_file;
        std::ofstream ofile;
        char obuff[300];
        uint j = 0;
        uint line_count = 0;

        // Syncfile parsing : current field and subfield in a line, and current (sub)field value is stored in temp
        uint field = 0;
        std::string temp = "";
        uint n_lines = 0;

        // Fields
        std::string contig = "";
        uint position = 0;
        char ref_allele = 'N';
        uint depth[2] = {0, 0};
        unsigned int pool[2][6] {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};

        // Base processing
        bool read_begin = false;
        bool next_indel = false;
        std::string tmp_indel_size = "";
        uint remaining_indel = 0;

        // Output
        std::string output_line = "";

        void add_counts_to_buffer(uint8_t pool_number);
        void process_line();
        void process_field();
        void process_base(char base, bool pool);
        void usage();
};
