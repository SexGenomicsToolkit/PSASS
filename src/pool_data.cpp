#include "pool_data.h"

void PoolBaseData::compute_total(){

    this->depth = std::accumulate(std::begin(this->nucleotides), std::end(this->nucleotides), uint16_t(0));
};


void PoolBaseData::compute_frequencies(){

    if (this->depth > 0) {

        for (auto i=0; i<6; ++i) this->frequencies[i] = float(this->nucleotides[i]) / float(this->depth);

    } else {

        for (auto i=0; i<6; ++i) this->frequencies[i] = 0;

    }
};


void PoolBaseData::compute_pi(){

    if (this->depth > 1) {

        this->pi = 1;
        for (auto i=0; i<6; ++i) this->pi -= this->frequencies[i] * this->frequencies[i];
        this->pi *= this->depth / (this->depth - 1);

    } else {

        this->pi = 0;

    }

};


void PoolBaseData::update(){

    this->compute_total();
    this->compute_frequencies();
    this->compute_pi();
};


std::string PoolBaseData::output_nucleotides(){

    std::stringstream output;

    for (auto i=0; i<6; ++i) {

        output << this->nucleotides[i];
        if (i < 5) output << "\t";

    }

    return output.str();
};


std::string PoolBaseData::output_frequencies(){

    std::stringstream output;
    output << std::fixed << std::setprecision(2);

    for (auto i=0; i<6; ++i) {

        output << this->frequencies[i];
        if (i < 5) output << "\t";

    }

    return output.str();
};
