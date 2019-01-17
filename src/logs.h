#pragma once
#include <fstream>
#include "utils.h"

class Logs {

    public:

        std::ofstream file;

        Logs() {}

        Logs(const std::string& output_prefix) {

            std::string file_path = output_prefix + ".log";
            this->file.open(file_path);

            if (not this->file.is_open()) {

                std::cerr << "Error: cannot open log file (" << file_path << ")." << std::endl;
                exit(1);

            }

        }

        template<typename T>
        void write(T line) {

            char logtime[DTTMSZ];
            this->file << "[" << print_time(logtime) << "]" << "  ";
            this->file << std::boolalpha << line;
            this->file << std::endl;

        }
};
