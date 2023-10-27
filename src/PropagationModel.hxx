#pragma once
#include <armadillo>
#include <functional>
#include <vector>
#include "WeightedEdgeGraph.h"

class PropagationModel{
    protected:
        std::function<double(double)> scaleFunction;
    public:
        //PropagationModel();
        //PropagationModel(std::function<double(double)> scaleFunction);
        virtual ~PropagationModel(){}
        //using the scale function as a parameter itself, dependency injection
        virtual arma::Col<double> propagate(arma::Col<double> input, double time) = 0;
        virtual arma::Col<double> propagationTerm(arma::Col<double> input, double time) = 0;

        //getters and setters
        std::function<double(double)> getScaleFunction(){return this->scaleFunction;}
        void setScaleFunction(std::function<double(double)> scaleFunction){this->scaleFunction = scaleFunction;}
};