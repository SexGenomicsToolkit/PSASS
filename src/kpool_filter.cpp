#include "kpool_filter.h"

void kpool_filter(Parameters& parameters) {

    std::ifstream input_file;
    input_file.open(parameters.input_file_path);

    std::string line;
    std::vector<std::string> fields;
    uint freq[2] {0, 0};
    uint lines = 0;

    std::getline(input_file, line);
    std::vector<std::string> header = split(line, "\t");

    std::string pool1_output_file_path = parameters.output_prefix + "_" + header[1] + ".tsv";
    std::ofstream pool1_output_file;

    std::string pool2_output_file_path = parameters.output_prefix + "_" + header[2] + ".tsv";
    std::ofstream pool2_output_file;

    bool output_pool1, output_pool2;

    parameters.pool_to_filter == "" or parameters.pool_to_filter == header[1] ? output_pool1 = true : output_pool1 = false;
    parameters.pool_to_filter == "" or parameters.pool_to_filter == header[2] ? output_pool2 = true : output_pool2 = false;

    if (output_pool1) {
        pool1_output_file.open(pool1_output_file_path);
        pool1_output_file << "sequence\t" << parameters.pool1_id << "\t" << parameters.pool2_id << "\n";
    }

    if (parameters.pool_to_filter == "" or parameters.pool_to_filter == header[2]) {
        pool2_output_file.open(pool2_output_file_path);
        pool2_output_file << "sequence\t" << parameters.pool1_id << "\t" << parameters.pool2_id << "\n";
    }

    while (std::getline(input_file, line)) {

        fields = split(line, "\t");
        freq[0] = std::stoi(fields[1]);
        freq[1] = std::stoi(fields[2]);

        if (output_pool1 and
            freq[0] >= parameters.min_kmer_presence_depth and
            freq[0] <= parameters.max_kmer_presence_depth and
            freq[1] <= parameters.max_kmer_absence_depth) {

            pool1_output_file << line << "\n";
        }

        if (output_pool2 and
            freq[1] >= parameters.min_kmer_presence_depth and
            freq[1] <= parameters.max_kmer_presence_depth and
            freq[0] <= parameters.max_kmer_absence_depth) {

            pool2_output_file << line << "\n";
        }

        if (lines % 25000000 == 0 and lines != 0) log("Processed " + std::to_string(lines / 1000000) + " M. lines.");
        ++lines;
    }
}
