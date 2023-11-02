#include "PropagationModelOriginal.hxx"
#include <armadillo>
#include <iostream>

PropagationModelOriginal::PropagationModelOriginal(const WeightedEdgeGraph* graph){
    this->scaleFunction = [](double time)-> double{return 0.5;};

    //pseudoinverse initialization
    std::vector<double> normalizationFactors(graph->getNumNodes(),0);
    for (int i = 0; i < graph->getNumNodes(); i++) {
        for(int j = 0; j < graph->getNumNodes();j++){
            normalizationFactors[i] += std::abs(graph->getEdgeWeight(i,j)); 
        }
    }
    arma::Mat<double> WtransArma = graph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
    
    arma::Mat<double> IdentityArma = arma::eye(graph->getNumNodes(),graph->getNumNodes());
    arma::Mat<double> temp = IdentityArma - WtransArma;
    
    //control determinant and invertibility, print warning if not invertible
    if(arma::det(temp) == 0){
        std::cout << "WARNING: The graph is not invertible, the pseudoinverse could lead to faulty results" << std::endl;
    }
    pseudoinverse = arma::pinv(temp);
}

PropagationModelOriginal::~PropagationModelOriginal(){
}

PropagationModelOriginal::PropagationModelOriginal(const WeightedEdgeGraph* graph,std::function<double(double)> scaleFunc):scaleFunction(scaleFunc){
    //pseudoinverse initialization
    std::vector<double> normalizationFactors(graph->getNumNodes(),0);
    for (int i = 0; i < graph->getNumNodes(); i++) {
        for(int j = 0; j < graph->getNumNodes();j++){
            normalizationFactors[i] += std::abs(graph->getEdgeWeight(i,j)); 
        }
    }
    arma::Mat<double> WtransArma = graph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
    
    arma::Mat<double> IdentityArma = arma::eye(graph->getNumNodes(),graph->getNumNodes());
    
    pseudoinverse = arma::pinv(IdentityArma - WtransArma);
}

arma::Col<double> PropagationModelOriginal::propagate(arma::Col<double> input, double time){
    return ( pseudoinverse * input * this->scaleFunction(time));
}

arma::Col<double> PropagationModelOriginal::propagationTerm(arma::Col<double> input, double time){
    //a propagation term doesn't exist in this case since it is a resolution of the system of equations
    return pseudoinverse * input * this->scaleFunction(time);
}
