#pragma once
#include <armadillo>
#include "DissipationModel.h"

class DissipationModelPow : public DissipationModel
{
    public:
        DissipationModelPow();
        DissipationModelPow(double power);
        ~DissipationModelPow();
        arma::Col<double> dissipate(arma::Col<double> input);
};