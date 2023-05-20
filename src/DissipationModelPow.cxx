#pragma once
#include <armadillo>
#include "DissipationModelPow.h"

DissipationModelPow::DissipationModelPow(){
    this->power = 2;
}

DissipationModelPow::DissipationModelPow(double power){
    this->power = power;
}

DissipationModelPow::~DissipationModelPow(){
}

arma::Col<double> DissipationModelPow::dissipate(arma::Col<double> input, double time){
    return input - pow(input,this->power);
}
