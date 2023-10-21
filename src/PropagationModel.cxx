#include "PropagationModel.hxx"
#include "armaUtilities.h"


PropagationModel::PropagationModel(){
    this->scaleFunction = [](double time)-> double{return 1;};
}

PropagationModel::PropagationModel(std::function<double(double)> scaleFunction){
    this->scaleFunction = scaleFunction;
}

PropagationModel::~PropagationModel(){}

arma::Col<double> PropagationModel::propagate(arma::Col<double> input, const WeightedEdgeGraph& graph, double time){
    
}

arma::Col<double> PropagationModel::propagationTerm(arma::Col<double> input, const WeightedEdgeGraph& graph, double time){
    
}


