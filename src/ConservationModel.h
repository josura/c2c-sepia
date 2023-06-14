#pragma once
#include <armadillo>
#include <functional>
#include "armaUtilities.h"

class ConservationModel{
    protected:
        std::function<double(double)> scaleFunction;
    public:
        ConservationModel();
        ConservationModel(std::function<double(double)> scaleFunction);
        virtual ~ConservationModel();
        //using the scale function as a parameter itself, dependency injection
        virtual arma::Col<double> conservate(arma::Col<double> input, arma::Col<double> inputDissipated,arma::Mat<double> Wstar, double time, std::vector<double> q = std::vector<double>());
        virtual arma::Col<double> conservationTerm(arma::Col<double> input,arma::Mat<double> Wstar, double time, std::vector<double> q = std::vector<double>());

        //getters and setters
        std::function<double(double)> getScaleFunction(){return this->scaleFunction;}
        void setScaleFunction(std::function<double(double)> scaleFunction){this->scaleFunction = scaleFunction;}
};