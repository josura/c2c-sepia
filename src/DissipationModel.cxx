#include "DissipationModel.h"

DissipationModel::DissipationModel(){
}

DissipationModel::~DissipationModel(){
}

arma::Col<double> DissipationModel::dissipatePow2Self(arma::Col<double> input){
    return input - pow(input,2);
}

