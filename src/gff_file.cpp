#include "gff_file.h"


GFFData::GFFData() {


}



bool GFFData::find_value(const std::string& field, const std::string& value) {

    if (this->field_values[field].find(value) != this->field_values[field].end()) return true;
    return false;
}



void GFFData::read_gff_file(std::ifstream& input_file, Logs& logs) {

    logs.write("Reading GFF file started.");

    std::string line = "", name = "", id = "", product = "", parent = "";
    std::vector<std::string> fields, infos;

    while (std::getline(input_file, line)) {

        if (not (line[0] == '#') and line.size() > 1) {

            name = "";
            id = "";
            product = "";
            parent = "";

            fields = split(line, "\t");

            if (fields.size() == 9) {  // Expected number of fields in GFF line

                infos = split(fields[8], ";");  // Last field contains non-standardized information

                if (find_value("gene", fields[2])) {  // Look for "gene" in type field

                    for (auto i: infos) {

                        if (find_value("name", i.substr(0, 4))) name = split(i, "=")[1];  // NCBI GFF have gene name in "Name" subfield of the last field
                        if (find_value("id", i.substr(0, 2))) id = split(i, "=")[1];  // Both NCBI and Non-NCBI GFF have an "ID" subfield which is also used to link other features

                    }

                    if (id != "") {

                        this->genes[id].contig = fields[0];
                        this->genes[id].start = fields[3];
                        this->genes[id].end = fields[4];
                        this->genes[id].name = name;
                        this->genes[id].id = id;

                        this->data[fields[0]].push_back(fields);  // Add the fields vector to GFF data, to create the GFF region later (used to check if current base is in gene / coding seq)

                    }


                } else if (find_value("transcript", fields[2])) {

                    for (auto i: infos) {

                        if (find_value("id", i.substr(0, 2))) id = split(i, "=")[1];  // Both NCBI and Non-NCBI GFF have an "ID" subfield which is also used to link other features
                        if (find_value("parent", i.substr(0, 6))) parent = split(i, "=")[1];  // Both NCBI and Non-NCBI GFF have a "Parent" subfield linking to gene

                    }

                    if (id != "" and parent != "") {

                        this->transcripts[id] = parent;

                    }

                }

                else if (find_value("CDS", fields[2])) {

                    for (auto i: infos) {

                        if (find_value("parent", i.substr(0, 6))) parent = split(i, "=")[1];
                        if (find_value("product", i.substr(0, 7))) product = split(i, "=")[1];

                    }

                    if (this->transcripts.find(parent) != this->transcripts.end()) {

                        this->genes[this->transcripts[parent]].product = product;
                        this->genes[this->transcripts[parent]].coding_length += std::stoul(fields[4]) - std::stoul(fields[3]) + 1;

                        this->data[fields[0]].push_back(fields);

                    }
                }
            }
        }
    }

    // Compute non-coding length for all genes
    for (auto& gene: this->genes) gene.second.noncoding_length = uint(std::stoul(gene.second.end) - std::stoul(gene.second.start) + 1 - gene.second.coding_length);

    logs.write("Reading GFF file ended without errors.");
    logs.write("Loaded <" + std::to_string(this->genes.size()) + "> genes from GFF.");
}



// Update the current contig per-base features. This is done for each contig to reduce memory overall.
void GFFData::new_contig(InputData& input_data, Logs& logs) {

    std::string feature = "", gene = "", parent = "";
    bool coding = false;
    std::vector<std::string> feature_info;

    this->contig.clear();

    for (auto fields: this->data[input_data.contig]) {

        feature = "";
        parent = "";
        feature_info = split(fields[8], ";");

        for (auto i: feature_info) {

            if (i.substr(0, 2) == "ID") feature = split(i, "=")[1];
            if (find_value("parent", i.substr(0, 6))) parent = split(i, "=")[1];

        }


        if (find_value("gene", fields[2])) {

            gene = feature;
            coding = false;

        } else if (find_value("CDS", fields[2])) {

            gene = this->transcripts[parent];
            coding = true;

        }

        for (auto i=std::stoul(fields[3]); i < std::stoul(fields[4]) + 1; ++i) this->contig[uint(i)] = std::pair<std::string, bool>(gene, coding);
    }

    logs.write("Loaded <" + std::to_string(this->contig.size()) + "> genic bases from GFF file for contig <" + input_data.contig + ">.");
}
