#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <map>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <sstream>

#define DTTMFMT "%Y-%m-%d %H:%M:%S"
#define DTTMSZ 21

// Output current date and time in format specified with DMTTMFMT and DTTMSZ
char* print_time (char *buff);

// Splits a std::string into a std::vector of std::strings according to a specified delimiter (default: \t)
std::vector<std::string> split(std::string str, const std::string delimiter);
