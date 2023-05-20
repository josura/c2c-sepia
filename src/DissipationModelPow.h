#pragma once
#include <armadillo>
#include "DissipationModel.h"

class DissipationModelPow : public DissipationModel
{
    private:
        double power;
    public:
        DissipationModelPow();
        DissipationModelPow(double power);
        ~DissipationModelPow();
        arma::Col<double> dissipate(arma::Col<double> input, double time);
};