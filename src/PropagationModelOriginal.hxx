#pragma once
#include <armadillo>
#include "PropagationModel.hxx"

class PropagationModelOriginal : public PropagationModel
{
    private:
        std::function<double(double)> scaleFunction;
        arma::dmat pseudoinverse;
    public:
        PropagationModelOriginal(const WeightedEdgeGraph* graph);
        PropagationModelOriginal(const WeightedEdgeGraph* graph,std::function<double(double)> scaleFunc);
        ~PropagationModelOriginal()override;
        arma::Col<double> propagate(arma::Col<double> input,double time)override;
        arma::Col<double> propagationTerm(arma::Col<double> input, double time)override;
        double getScale(double time){return scaleFunction(time);}
};