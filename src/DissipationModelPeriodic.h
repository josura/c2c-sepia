#pragma once
#include <armadillo>
#include "DissipationModel.h"

class DissipationModelPeriodic : public DissipationModel
{
    public:
        DissipationModelPeriodic();
        DissipationModelPeriodic(double phase, double period, double amplitude);
        ~DissipationModelPeriodic();
        arma::Col<double> dissipate(arma::Col<double> input);
};