#pragma once
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <vector>
#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include "parameters.h"
#include "utils.h"

// Simple structure holding all information about an input file
struct inputFile {
    htsFile *sam;  // Main file descriptor (for the alignment file)
    hts_idx_t *idx;  // Index file descriptor
    sam_hdr_t *header;  // Header information read directly from main file
    uint16_t file_n;  // Input file number
};

// Open an alignment file in a format-agnostic way and fill an inputFile object with all the information
int open_input(std::string &fn_in, inputFile *file, std::string &reference, uint16_t file_n);

// Compute nucleotide depths for a region in a single alignment file
int process_file(inputFile* input, char *region, std::vector<std::vector<uint16_t>>& depths, uint min_qual=0);

// Main pileup function
int pileup(Parameters& parameters);
