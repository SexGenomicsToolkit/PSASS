#pragma once
#include <set>
#include <string>
#include <unordered_map>
#include "input_data.h"
#include "utils.h"


struct Gene {

    std::string contig = "";
    std::string start = "";
    std::string end = "";
    std::string name = "";
    std::string id = "";
    std::string product = "";
    uint coding_length = 0;
    uint noncoding_length = 0;
    uint depth[6] {0, 0, 0, 0, 0, 0}; // Coding male, Non-coding male, Coding female, Non-coding female, Total male, Total female
    uint snps[6] {0, 0, 0, 0, 0, 0}; // Coding male, Non-coding male, Coding female, Non-coding female, Total male, Total female
};


class GFFData {

    public:

        std::unordered_map<std::string, std::set<std::string>> field_values= {{"gene",{"gene"}},
                                                                              {"transcript",{"mRNA", "transcript"}},
                                                                              {"CDS",{"CDS"}},
                                                                              {"name",{"Name"}},
                                                                              {"id",{"ID"}},
                                                                              {"parent",{"Parent"}},
                                                                              {"product",{"product"}}};

        std::unordered_map<std::string, std::vector<std::vector<std::string>>> data;
        std::unordered_map<std::string, Gene> genes;
        std::unordered_map<uint, std::pair<std::string, bool>> contig;
        std::unordered_map<std::string, std::string> transcripts;

        GFFData();
        void read_gff_file(std::istream& input_file);
        void new_contig(InputData& input_data);

    private:

        bool find_value(const std::string& field, const std::string& value);

};


