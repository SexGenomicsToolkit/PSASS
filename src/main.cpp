#include "arg_parser.h"
#include "analysis.h"
#include "psass.h"
#include <chrono>


int main(int argc, char *argv[]) {

//    std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();

    Psass psass(argc, argv);
    psass.run();

    std::cout << "Done" << std::endl;


//    std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
//    write_log("Analysis completed. Processed ", parameters.log_file, true, false);
//    write_log(n_lines, parameters.log_file, false, false);
//    write_log(" lines in ", parameters.log_file, false, false);
//    write_log(std::chrono::duration_cast<std::chrono::seconds>(t_end - t_begin).count(), parameters.log_file, false, false);
//    write_log(" seconds.", parameters.log_file, false, false);

    return 0;
}
