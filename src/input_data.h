#pragma once
#include <string>


class InputData {

    public:

        // Efficient file reading parameters
        char buff[2048];
        int buff_size = sizeof(this->buff);
        long k = 0;

        // Syncfile parsing : current field and subfield in a line, and current (sub)field value is stored in temp
        uint field = 0;
        uint subfield = 0;
        std::string temp = "";

        // Contig and position information
        std::string contig = "";
        std::string current_contig = "";
        uint position = 0;

        InputData() {}
};

