#include "DissipationModelPow.h"

DissipationModelPow::DissipationModelPow(){
    this->power = 2;
    this->numEl = 0;
}

DissipationModelPow::DissipationModelPow(double power){
    this->power = power;
    this->numEl = 0;
}

DissipationModelPow::~DissipationModelPow(){
}

arma::Col<double> DissipationModelPow::dissipate(arma::Col<double> input, double time){
    return input - pow(input,this->power);
}

arma::Col<double> DissipationModelPow::dissipationTerm(arma::Col<double> input, double time){
    return pow(input,this->power);
}
