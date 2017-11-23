#include "arg_parser.h"
#include "analysis.h"


int main(int argc, char *argv[]) {

    ArgParser cmd_options(argc, argv);

    Parameters parameters;

    cmd_options.set_parameters(parameters);

    cmd_options.print_parameters();

    std::cout << "Analyzing" << std::endl;
    analysis(parameters);

    return 0;
}
