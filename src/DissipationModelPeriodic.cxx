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

DissipationModelPeriodic::DissipationModelPeriodic(arma::Col<double> phases, arma::Col<double> periods, arma::Col<double> amplitudes){
    this->numEl = phases.n_elem;
    this->phases = phases;
    this->periods = periods;
    this->amplitudes = amplitudes;
}

DissipationModelPeriodic::~DissipationModelPeriodic(){
}

arma::Col<double> DissipationModelPeriodic::dissipate(arma::Col<double> input, double time){
    arma::Col<double> output = arma::Col<double>(this->numEl);
    for(int i = 0; i < this->numEl; i++){
        output(i) = input(i) - this->amplitudes(i)*sin(2*arma::datum::pi/this->periods(i)*time + this->phases(i));
    }
    return output;
}