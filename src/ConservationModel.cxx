#include "ConservationModel.h"
#include "armaUtilities.h"


ConservationModel::ConservationModel(){
    this->scaleFunction = [](double time)-> double{return 0.5;};
}

ConservationModel::ConservationModel(std::function<double(double)> scaleFunction){
    this->scaleFunction = scaleFunction;
}

ConservationModel::~ConservationModel(){}

arma::Col<double> ConservationModel::conservate(arma::Col<double> input, arma::Col<double> inputDissipated, arma::Mat<double> Wstar,double time, std::vector<double> q){
    //if q is empty, then we assume that all the values in q are 1, that means all the weights of the edges are considered of the same importance and 
    // and the perturbation is completely passed down the network 
    if (q.size()) {
        if (q.size() == input.n_elem) {
            //convert q vector to arma vector
            arma::Col<double> qArma = vectorToArmaColumn(q);
            arma::Col<double> outputArma = inputDissipated -  scaleFunction(time) * Wstar * qArma * input;
            return outputArma;
        } else{
            throw std::invalid_argument("q vector is not of the same size as input vector. abort");
        }
    }
    else {
        //since in the case of vector values with all q values equal to 1
        arma::Col<double> qOnes = arma::ones<arma::Col<double>>(input.n_elem);
        arma::Col<double> outputArma =  inputDissipated -  scaleFunction(time)* Wstar * qOnes * input;
        return outputArma;
    }
}

arma::Col<double> ConservationModel::conservationTerm(arma::Col<double> input, arma::Mat<double> Wstar, double time, std::vector<double> q){
    if (q.size()) {
        if (q.size() == input.n_elem) {
            //convert q vector to arma vector
            arma::Col<double> qArma = vectorToArmaColumn(q);
            arma::Col<double> outputArma =  scaleFunction(time) * Wstar * qArma * input;
            return outputArma;
        } else{
            throw std::invalid_argument("q is not of the same size as input vector. abort");
        }
    } else{
        arma::Col<double> qOnes = arma::ones<arma::Col<double>>(input.n_elem);
        arma::Col<double> outputArma =  scaleFunction(time) * Wstar * qOnes * input;
        return outputArma;
    }
}


