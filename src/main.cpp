#include "arg_parser.h"
#include "analysis.h"


int main(int argc, char *argv[]) {

    ArgParser cmd_options(argc, argv);

    Parameters parameters;

    cmd_options.set_parameters(parameters);

    cmd_options.print_parameters();

    //    for (auto gene: genes) {
//        std::cout << gene.first << "," << gene.second.name << "," << gene.second.contig << "," << gene.second.start << "," << gene.second.end << "," << gene.second.product << ","
//                  << gene.second.coverage[0] << "," << gene.second.coverage[1] << "," << gene.second.coverage[2] << "," << gene.second.coverage[3] << ","
//                  << gene.second.snps[0] << "," << gene.second.snps[1] << "," << gene.second.snps[2] << "," << gene.second.snps[3] << "\n";
//    }

//    for (auto region: regions) {
//        for (auto pos: region.second) std::cout << region.first << "," << pos.first << "," << pos.second << "\n";
//    }

    std::cout << "Analyzing" << std::endl;
    analysis(parameters);

    return 0;
}
