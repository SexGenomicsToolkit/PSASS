#include "arg_parser.h"
#include "parameters.h"
#include "pileup_converter.h"
#include "pileup.h"
#include "psass.h"


int main(int argc, char *argv[]) {

    Parameters parameters = parse_args(argc, argv);

    if (parameters.command == "analyze") {

        Psass psass(parameters);
        psass.run();

    } else if (parameters.command == "convert") {

        PileupConverter pileup_converter(parameters);
        pileup_converter.run();

    } else if (parameters.command == "pileup") {

        pileup(parameters);

    }

    return 0;
}
