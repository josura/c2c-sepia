#include "DissipationModelScaled.h"

DissipationModelScaled::DissipationModelScaled(){
    this->scaleFunction = [](double time)-> double{return 0.5;};
    this->numEl = 0;
}

DissipationModelScaled::DissipationModelScaled(std::function<double(double)> scaleFun){
    this->scaleFunction = scaleFun;
    this->numEl = 0;
}

DissipationModelScaled::~DissipationModelScaled(){
}

arma::Col<double> DissipationModelScaled::dissipate(arma::Col<double> input, double time){
    return input - this->scaleFunction(time)*input;
}

arma::Col<double> DissipationModelScaled::dissipationTerm(arma::Col<double> input, double time){
    return this->scaleFunction(time)*input;
}
