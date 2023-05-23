#pragma once
#include <armadillo>
#include "DissipationModel.h"

class DissipationModelRandom : public DissipationModel
{
    private:
        double rangeMin;
        double rangeMax;
    public:
        DissipationModelRandom();
        DissipationModelRandom(double rangeMin, double rangeMax);
        ~DissipationModelRandom();
        arma::Col<double> dissipate(arma::Col<double> input, double time)override;
        arma::Col<double> dissipationTerm(arma::Col<double> input, double time)override;
};