

#include "Computation.h"
#include "Matrix.h"
#include <map>
#include <vector>


Computation::Computation(){
    input=std::vector<double>();
    output=std::vector<double>();
    W = Matrix<double>();
    Wstar = Matrix<double>();
    geneMapToNode = std::map<std::string,int>();
    cellTypes = std::vector<std::string>();
}

Computation::Computation(std::vector<double>& _input){
    input=_input;
    output=std::vector<double>();
    W = Matrix<double>();
    Wstar = Matrix<double>();
    geneMapToNode = std::map<std::string,int>();
    cellTypes = std::vector<std::string>();
}

Computation::Computation(std::vector<double>& _input, Matrix<double>& _W){
    input=_input;
    output=std::vector<double>();
    W = _W;
    Wstar = Matrix<double>();
    geneMapToNode = std::map<std::string,int>();
    cellTypes = std::vector<std::string>();
}

Computation::Computation(std::vector<double>& _input, Matrix<double>& _W,std::vector<std::string>& _cellTypes){
    input=_input;
    output=std::vector<double>();
    W = _W;
    Wstar = Matrix<double>();
    geneMapToNode = std::map<std::string,int>();
    cellTypes = _cellTypes;
}