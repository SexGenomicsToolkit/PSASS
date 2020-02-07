#include "kpool_merge.h"


void kpool_merge(Parameters& parameters) {

    // Open input files
    std::ifstream table_1_file, table_2_file;
    table_1_file.open(parameters.table1_file_path);
    table_2_file.open(parameters.table1_file_path);

    IndexFiles index_files;

    for (auto index: indexes) {
        index_files[index] = {parameters.tmp_file_prefix + index + ".1.kpool.tmp", parameters.tmp_file_prefix + index + ".2.kpool.tmp"};
    }

    log("Indexing kmer tables");
    std::mutex log_mutex;
    std::thread table_1_thread(index_file, std::ref(table_1_file), std::ref(index_files), std::ref(log_mutex), 1);
    std::thread table_2_thread(index_file, std::ref(table_2_file), std::ref(index_files), std::ref(log_mutex), 2);
    table_1_thread.join();
    table_2_thread.join();

    uint64_t kmer_count = 0;

    log("Merging indexed files");
    std::ofstream output_file;
    output_file.open(parameters.output_file_path);

    output_file << "sequence" << "\t" << parameters.pool1_id << "\t" << parameters.pool2_id << std::endl;

    for (auto index: index_files) {
        process_batch(index.second, output_file, kmer_count);
    }

   log("Merge ended successfully. Final kmer count : " + std::to_string(kmer_count));
}



void index_file(std::ifstream& input_file, IndexFiles& index_file_paths, std::mutex& log_mutex, uint8_t pool) {

    std::unordered_map<std::string, std::ofstream> index_files;
    for (auto index: index_file_paths) index_files[index.first].open(index.second[pool - 1]);

    std::string line;
    uint lines = 0;

    while (std::getline(input_file, line)) {

        if (lines % 25000000 == 0 and lines != 0) {
            log_mutex.lock();
            log("Indexed " + std::to_string(lines / 1000000) + " M. kmers in pool " + std::to_string(pool));
            log_mutex.unlock();
        }
        ++lines;
        index_files[line.substr(0, 3)] << line << "\n";
    }

    for (auto index: index_file_paths) index_files[index.first].close();
}



void process_batch(std::vector<std::string>& tmp_files, std::ofstream& output_file, uint64_t& kmer_count) {

    std::unordered_map<std::string, std::pair<uint32_t, uint32_t>> table;

    std::ifstream pool_1_file, pool_2_file;
    pool_1_file.open(tmp_files[0]);
    pool_2_file.open(tmp_files[1]);

    std::string line;
    std::vector<std::string> fields;

    log("Processing index file <" + tmp_files[0] + ">");

    while(std::getline(pool_1_file, line)) {
        fields = split(line, "\t");
        table[fields[0]].first = std::stoi(fields[1]);
    }

    log("Processing index file <" + tmp_files[1] + ">");

    while(std::getline(pool_2_file, line)) {
        fields = split(line, "\t");
        table[fields[0]].second = std::stoi(fields[1]);
    }

    kmer_count += table.size();

    for (auto kmer: table) {
        output_file << kmer.first << "\t" << kmer.second.first << "\t" << kmer.second.second << "\n";
    }

    log("Removing temporary index files");
    remove(tmp_files[0].c_str());
    remove(tmp_files[1].c_str());
}
