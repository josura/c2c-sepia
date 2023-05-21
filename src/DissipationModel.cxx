#include "DissipationModel.h"

DissipationModel::DissipationModel(){
}

DissipationModel::~DissipationModel(){
}


arma::Col<double> DissipationModel::dissipate(arma::Col<double> input,double time){
    return input;
}
arma::Col<double> DissipationModel::dissipationTerm(arma::Col<double> input, double time){
    return input;
}

arma::Col<double> DissipationModel::dissipatePow2Self(arma::Col<double> input){
    return input - pow(input,2);
}

arma::Col<double> DissipationModel::dissipateSelfScaled(arma::Col<double> input, double scale){
    return input - scale*input;
}

arma::Col<double> DissipationModel::dissipateSelfScaled(arma::Col<double> input, arma::Col<double> scales){
    return input - scales%input;
}

arma::Col<double> DissipationModel::dissipateSelfPeriodic(arma::Col<double> input, double period, double amplitude, double phase){
    return input;
}

arma::Col<double> DissipationModel::dissipateSelfPeriodic(arma::Col<double> input, arma::Col<double> periods, arma::Col<double> amplitudes, arma::Col<double> phases){
    return input;
}

arma::Col<double> DissipationModel::dissipateSelfRandom(arma::Col<double> input, double min, double max){
    return input;
}