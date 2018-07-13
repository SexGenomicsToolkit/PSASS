#include "gff_file.h"


void read_gff_file(std::ifstream& input_file, std::unordered_map<std::string, std::unordered_map<uint, std::string>>& regions, std::unordered_map<std::string, Gene>& genes) {

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

                    genes[gene].contig = fields[0];
                    genes[gene].start = fields[3];
                    genes[gene].end = fields[4];
                    genes[gene].name = gene;

                } else if (fields[2] == "exon") {

                    for (auto i: infos) {
                        if (i.substr(0, 4) == "gene") gene = split(i, "=")[1];
                        else if (i.substr(0, 4) == "product") product = split(i, "=")[1];
                    }

                    if (gene != "" and (genes.find(gene) != genes.end())) {
                        genes[gene].product = product;
                        for (int i=std::stoi(fields[3]); i < std::stoi(fields[4]) + 1; ++i) regions[fields[0]][i] = gene;
                    }
                }
            }
        }
    }
}

