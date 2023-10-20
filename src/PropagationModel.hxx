#pragma once
#include <armadillo>
#include <functional>

class PropagationModel{
    protected:
        std::function<double(double)> scaleFunction;
    public:
        PropagationModel();
        PropagationModel(std::function<double(double)> scaleFunction);
        virtual ~PropagationModel();
        //using the scale function as a parameter itself, dependency injection
        virtual arma::Col<double> propagate(arma::Col<double> input, arma::Col<double> inputDissipated,arma::Mat<double> Wstar, double time, std::vector<double> q = std::vector<double>());
        virtual arma::Col<double> propagationTerm(arma::Col<double> input,arma::Mat<double> Wstar, double time, std::vector<double> q = std::vector<double>());

        //getters and setters
        std::function<double(double)> getScaleFunction(){return this->scaleFunction;}
        void setScaleFunction(std::function<double(double)> scaleFunction){this->scaleFunction = scaleFunction;}
};