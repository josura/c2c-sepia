#include "DissipationModelScaled.h"
#include <iostream>

DissipationModelScaled::DissipationModelScaled(){
    this->scaleFunction = [](double time)-> double{return 0.5;};
    this->numEl = 0;
}

DissipationModelScaled::~DissipationModelScaled(){
}

DissipationModelScaled::DissipationModelScaled(std::function<double(double)> scaleFun){
    this->scaleFunction = scaleFun;
    this->numEl = 0;
    std::cout << "DissipationModelScaled constructor" << std::endl;
}

arma::Col<double> DissipationModelScaled::dissipate(arma::Col<double> input, double time){
    return input - (this->scaleFunction(time)*input);
    std::cout << "dissipate from scaled" << std::endl;
}

arma::Col<double> DissipationModelScaled::dissipationTerm(arma::Col<double> input, double time){
    return this->scaleFunction(time)*input;
    std::cout << "dissipate term from scaled" << std::endl;
}
