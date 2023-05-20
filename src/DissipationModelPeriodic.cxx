#include "DissipationModelPeriodic.h"

DissipationModelPeriodic::DissipationModelPeriodic(){
    numEl = 1;
    this->phases = arma::Col<double>(numEl);
    this->periods = arma::Col<double>(numEl);
    this->amplitudes = arma::Col<double>(numEl);
}

DissipationModelPeriodic::DissipationModelPeriodic(int numEl,double phase, double period, double amplitude){
    this->numEl = numEl;
    this->phases = arma::Col<double>(numEl,arma::fill::scalar_holder(phase));
    this->periods = arma::Col<double>(numEl,arma::fill::scalar_holder(period));
    this->amplitudes = arma::Col<double>(numEl,arma::fill::scalar_holder(amplitude));
}