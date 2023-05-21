#include "DissipationModelRandom.h"

DissipationModelRandom::DissipationModelRandom(){
    this->rangeMin = 0;
    this->rangeMax = 1;
    this->numEl = 1;
}

DissipationModelRandom::DissipationModelRandom(int numEl,double rangeMin, double rangeMax){
    this->rangeMin = rangeMin;
    this->rangeMax = rangeMax;
    this->numEl = numEl;
}

DissipationModelRandom::~DissipationModelRandom(){
}

arma::Col<double> DissipationModelRandom::dissipate(arma::Col<double> input, double time){
    arma::Col<double> output = arma::Col<double>(this->numEl);
    for(int i = 0; i < this->numEl; i++){
        output(i) = input(i) - input(i)*(this->rangeMin + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(this->rangeMax-this->rangeMin))));
    }
    return output;
}

arma::Col<double> DissipationModelRandom::dissipateTerm(arma::Col<double> input, double time){
    arma::Col<double> output = arma::Col<double>(this->numEl);
    for(int i = 0; i < this->numEl; i++){
        output(i) = input(i)*(this->rangeMin + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(this->rangeMax-this->rangeMin))));
    }
    return output;
}