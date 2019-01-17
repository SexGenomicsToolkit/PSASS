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

    time_t t = time(nullptr);
    strftime(buff, DTTMSZ, DTTMFMT, localtime (&t));
    return buff;
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
inline int fast_stoi(const char* str) {

    int val = 0;
    while( *str ) {
        val = val*10 + (*str++ - '0');
    }
    return val;
}
