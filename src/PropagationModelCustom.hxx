#pragma once
#include <armadillo>
#include "PropagationModel.hxx"

class PropagationModelCustom : public PropagationModel
{
    private:
        std::function<double(double)> scaleFunction;
        arma::dmat Wmat;
    public:
        PropagationModelCustom(const WeightedEdgeGraph* graph);
        PropagationModelCustom(const WeightedEdgeGraph* graph,std::function<double(double)> scaleFunc);
        ~PropagationModelCustom()override;
        arma::Col<double> propagate(arma::Col<double> input,double time)override; // add additional parameters, but remember to change the main accordingly
        arma::Col<double> propagationTerm(arma::Col<double> input, double time)override;
        double getScale(double time){return scaleFunction(time);}
};