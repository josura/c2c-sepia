#pragma once
#include <armadillo>
#include "DissipationModel.h"

class DissipationModelRandom : public DissipationModel
{
    private:
        double rangeMin;
        double rangeMax;
        int numEl;
    public:
        DissipationModelRandom();
        DissipationModelRandom(int numEl,double rangeMin, double rangeMax);
        ~DissipationModelRandom();
        arma::Col<double> dissipate(arma::Col<double> input, double time);
};