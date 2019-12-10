#include "psass.h"
#include "pileup_converter.h"


int main(int argc, char *argv[]) {

    if (argc < 2 or (std::string(argv[1]) != "analyze" and std::string(argv[1]) != "convert")) {

        std::cerr << std::endl << "Usage: poolsex <command> [options]" << std::endl << std::endl;
        std::cerr << "Commands:" << std::endl;
        std::cerr << "        analyze    Compute metrics from either the resulting file from \"convert\" or a sync file from popoolation2" << std::endl;
        std::cerr << "        convert    Convert a pileup file to a synchronized pool file. Accepts input from stdin with \"-\"" << std::endl;
        std::cerr << std::endl;

    } else if (std::string(argv[1]) == "analyze") {

        Psass psass(argc, argv);
        psass.run();

    } else if (std::string(argv[1]) == "convert") {

        PileupConverter pileup_converter(argc, argv);
        pileup_converter.run();

    }

    return 0;
}
