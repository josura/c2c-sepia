#include "DissipationModel.h"

DissipationModel::DissipationModel(){
}

DissipationModel::~DissipationModel(){
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
