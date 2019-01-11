#include "pair_data.h"

PairBaseData::PairBaseData() {

}


void PairBaseData::compute_average_freq() {

    for (auto i=0; i<6; ++i) this->average_frequencies[i] = (this->pool1.frequencies[i] + this->pool2.frequencies[i]) / 2;
};


void PairBaseData::compute_total_pi() {

    this->total_pi = 1;
    for (auto i=0; i<6; ++i) this->total_pi -= this->average_frequencies[i] * this->average_frequencies[i];
    float min_total = std::min(this->pool1.depth, this->pool2.depth);
    this->total_pi *= min_total / (min_total - 1);
};


void PairBaseData::compute_within_pi() {

    this->within_pi = (this->pool1.pi + this->pool2.pi) / 2;
}


void PairBaseData::compute_fst() {

    if (this->total_pi > 0) {
        this->fst = (this->total_pi - this->within_pi) / this->total_pi;
    } else {
        this->fst = 0;
    }
}


void PairBaseData::update() {

    this->pool1.update();
    this->pool2.update();
    this->compute_average_freq();
    this->compute_total_pi();
    this->compute_within_pi();
    this->compute_fst();
};
