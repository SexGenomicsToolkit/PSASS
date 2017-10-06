#include "output.h"


void output_data(variants &data, std::ofstream& file) {

//    file << "Contig" << "\t" << "Position" << "\t" << "M_cov" << "\t" << "F_cov" << "\t" << "M_h" << "\t" << "F_h" << "\n";
//    for (auto contig : data) {
//        for (uint i=0; i<contig.second.size(); ++i) {
//            file << contig.first << "\t" << i+1 << "\t" << contig.second[i].male_coverage << "\t" << contig.second[i].female_coverage << "\t";
//            file << contig.second[i].male_heterozygote << "\t" << contig.second[i].female_heterozygote << "\n";
//        }
//    }

    file << "Contig" << "\t" << "Position" << "\t"
         << "M_A" << "\t" << "M_T" << "\t" << "M_G" << "\t" << "M_C" << "\t" << "M_I" << "\t"
         << "F_A" << "\t" << "F_T" << "\t" << "F_G" << "\t" << "F_C" << "\t" << "F_I" << "\n";

    for (auto contig : data) {
        for (uint i=0; i<contig.second.size(); ++i) {
            file << contig.first << "\t" << i+1 << "\t"
                 << contig.second[i].male_bases[0] << "\t"
                 << contig.second[i].male_bases[1] << "\t"
                 << contig.second[i].male_bases[2] << "\t"
                 << contig.second[i].male_bases[3] << "\t"
                 << contig.second[i].male_bases[4] << "\t"
                 << contig.second[i].female_bases[0] << "\t"
                 << contig.second[i].female_bases[1] << "\t"
                 << contig.second[i].female_bases[2] << "\t"
                 << contig.second[i].female_bases[3] << "\t"
                 << contig.second[i].female_bases[4] << "\n";
        }
    }
}
