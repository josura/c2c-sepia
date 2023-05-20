#pragma once
#include <armadillo>
#include "DissipationModel.h"

class DissipationModelPow : public DissipationModel
{
    DissipationModelPow();
    DissipationModelPow(double power);
    ~DissipationModelPow();
    public:
        arma::Col<double> dissipate(arma::Col<double> input);
};