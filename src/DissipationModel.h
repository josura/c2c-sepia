#pragma once
#include <armadillo>

class DissipationModel{
    DissipationModel();
    ~DissipationModel();
    public:
        arma::Col<double> dissipatePow2Self(arma::Col<double> input);
};