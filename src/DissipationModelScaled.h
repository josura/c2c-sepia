#pragma once
#include <armadillo>
#include "DissipationModel.h"

class DissipationModelScaled : public DissipationModel
{
    private:
        std::function<double(double)> scaleFunction;
    public:
        // default constructor uses scaleFunction = 0.5
        DissipationModelScaled();
        DissipationModelScaled(std::function<double(double)> scaleFunc);
        ~DissipationModelScaled()override;
        arma::Col<double> dissipate(arma::Col<double> input, double time)override;
        arma::Col<double> dissipationTerm(arma::Col<double> input, double time)override;
        double getScale(double time){return scaleFunction(time);}
};