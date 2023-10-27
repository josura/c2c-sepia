#pragma once
#include <armadillo>
#include "PropagationModel.hxx"

class PropagationModelNeighbors : public PropagationModel
{
    private:
        std::function<double(double)> scaleFunction;
        arma::dmat Wmat;
    public:
        PropagationModelNeighbors(const WeightedEdgeGraph* graph);
        PropagationModelNeighbors(const WeightedEdgeGraph* graph,std::function<double(double)> scaleFunc);
        ~PropagationModelNeighbors()override;
        arma::Col<double> propagate(arma::Col<double> input,double time)override;
        arma::Col<double> propagationTerm(arma::Col<double> input, double time)override;
        double getScale(double time){return scaleFunction(time);}
};