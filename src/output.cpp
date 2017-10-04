#include "output.h"


void output_data(variants &data, std::ofstream& file) {

    file << "Contig" << "\t" << "Position" << "\t" << "M_cov" << "\t" << "F_cov" << "\t" << "M_h" << "\t" << "F_h" << "\n";
    for (auto contig : data) {
        for (uint i=0; i<contig.second.size(); ++i) {
            file << contig.first << "\t" << i+1 << "\t" << contig.second[i].male_coverage << "\t" << contig.second[i].female_coverage << "\t";
            file << contig.second[i].male_heterozygote << "\t" << contig.second[i].female_heterozygote << "\n";
        }
    }
}
