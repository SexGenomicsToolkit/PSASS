#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <map>

#define DTTMFMT "%Y-%m-%d %H:%M:%S"
#define DTTMSZ 21


// Output current date and time in format specified in utils.h
inline char* print_time (char *buff) {

    time_t t = time(0);
    strftime(buff, DTTMSZ, DTTMFMT, localtime (&t));
    return buff;
}



// Output a message in log file
template<typename T>
inline void write_log(T line, std::ofstream& log_file, bool timestamp, bool newline) {

    char logtime[DTTMSZ];
    if (timestamp) log_file << "[" << print_time(logtime) << "]" << "    ";
    log_file << std::boolalpha << line;
    if (newline) log_file << std::endl;
}



// Splits a std::string into a std::vector of std::strings according to a specified delimiter (default: \t)
inline std::vector<std::string> split(std::string str, const std::string delimiter) {

    std::vector<std::string> output;
    size_t pos;

    while ((pos = str.find(delimiter)) != std::string::npos){

        output.push_back(str.substr(0, pos));
        str.erase(0, pos + delimiter.length());
    }

    output.push_back(str.substr(0, pos));

    return output;
}



// Faster string to int conversion
inline uint64_t fast_stoi(const char* str) {

    uint64_t val = 0;
    while( *str ) {
        val = val*10 + (*str++ - '0');
    }
    return val;
}
