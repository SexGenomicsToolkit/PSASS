#pragma once
#include <mutex>
#include <thread>
#include <unordered_map>
#include "parameters.h"
#include "utils.h"

const std::vector<std::string> indexes {"AAA", "AAT", "AAG", "AAC",
                                        "ATA", "ATT", "ATG", "ATC",
                                        "ACA", "ACT", "ACG", "ACC",
                                        "AGA", "AGT", "AGG", "AGC",
                                        "TAA", "TAT", "TAG", "TAC",
                                        "TTA", "TTT", "TTG", "TTC",
                                        "TGA", "TGT", "TGG", "TGC",
                                        "TCA", "TCT", "TCG", "TCC",
                                        "CAA", "CAT", "CAG", "CAC",
                                        "CTA", "CTT", "CTG", "CTC",
                                        "CGA", "CGT", "CGG", "CGC",
                                        "CCA", "CCT", "CCG", "CCC",
                                        "GAA", "GAT", "GAG", "GAC",
                                        "GTA", "GTT", "GTG", "GTC",
                                        "GGA", "GGT", "GGG", "GGC",
                                        "GCA", "GCT", "GCG", "GCC",
                                       };

typedef std::unordered_map<std::string, std::vector<std::string>> IndexFiles;

void kpool_merge(Parameters& parameters);

void index_file(std::ifstream& input_file, IndexFiles& index_file_paths, std::mutex& log_mutex, uint8_t suffix);

void process_batch(std::vector<std::string>& tmp_files, std::ofstream& output_file, uint64_t& kmer_count);
