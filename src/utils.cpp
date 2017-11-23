#include "utils.h"

uint32_t fast_stoi(const char* str){

    uint32_t val = 0;
    while( *str ) {
        val = val*10 + (*str++ - '0');
    }
    return val;
}
