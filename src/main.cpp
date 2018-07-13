#include "arg_parser.h"
#include "analysis.h"
#include <chrono>


int main(int argc, char *argv[]) {

    std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();

    ArgParser cmd_options(argc, argv);

    Parameters parameters;

    cmd_options.set_parameters(parameters);

    cmd_options.print_parameters();

    std::cout << "Analyzing" << std::endl;
    uint n_lines = analysis(parameters);

    std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
    write_log("Analysis completed. Processed ", parameters.log_file, true, false);
    write_log(n_lines, parameters.log_file, false, false);
    write_log(" lines in ", parameters.log_file, false, false);
    write_log(std::chrono::duration_cast<std::chrono::seconds>(t_end - t_begin).count(), parameters.log_file, false, false);
    write_log(" seconds.", parameters.log_file, false, false);

    return 0;
}
