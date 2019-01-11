#include "gff_file.h"


void read_gff_file(std::ifstream& input_file, std::unordered_map<std::string, std::vector<std::vector<std::string>>>& file, std::unordered_map<std::string, Gene>& genes) {

    std::string line, gene, product;
    std::vector<std::string> fields, infos;

    while (std::getline(input_file, line)) {

        if (not (line[0] == '#') and line.size() > 1) {

            gene = "";
            product = "";

            fields = split(line, "\t");

            if (fields.size() == 9) {

                infos = split(fields[8], ";");

                if (fields[2] == "gene") {

                    for (auto i: infos) {
                        if (i.substr(0, 4) == "Name") gene = split(i, "=")[1];
                    }

                    if (gene == "") {
                        for (auto i: infos) {
                            if (i.substr(0, 2) == "ID") gene = split(i, "=")[1];
                        }
                    }

                    if (gene != "") {
                        genes[gene].contig = fields[0];
                        genes[gene].start = fields[3];
                        genes[gene].end = fields[4];
                        genes[gene].name = gene;
                        file[fields[0]].push_back(fields);
                    }


                } else if (fields[2] == "exon") {

                    for (auto i: infos) {
                        if (i.substr(0, 4) == "gene") gene = split(i, "=")[1];
                        else if (i.substr(0, 7) == "product") product = split(i, "=")[1];
                    }

                    if (gene == "") {
                        for (auto i: infos) {
                            if (i.substr(0, 2) == "ID") gene = split(i, "=")[1];
                        }
                    }

                    if (gene != "" and (genes.find(gene) != genes.end())) {
                        genes[gene].product = product;
                        genes[gene].coding_length += std::stoul(fields[4]) - std::stoul(fields[3]);
//                        file[fields[0]].push_back(fields);
                    }
                }
            }
        }
    }
}

