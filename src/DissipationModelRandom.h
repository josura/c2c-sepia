#pragma once
#include <armadillo>
#include "DissipationModel.h"

class DissipationModelRandom : public DissipationModel
{
    public:
        DissipationModelRandom();
        DissipationModelRandom(double rangeMin, double rangeMax);
        ~DissipationModelRandom();
        arma::Col<double> dissipate(arma::Col<double> input);
};