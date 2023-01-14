#pragma once

#include "Matrix.h"
#include "WeightedEdgeGraph.h"
#include <map>
#include <string>
#include <tuple>
#include <vector>

class Computation{
    private:
        std::vector<double> input,output;
        WeightedEdgeGraph* metapathway;
        WeightedEdgeGraph* augmentedMetapathway;
        std::map<std::string, int> geneMapToNode;
        std::vector<std::string> cellTypes;
    public:
        Computation();
        Computation(const std::vector<double>& _input);
        Computation(const std::vector<double>& _input, const Matrix<double>& _W);
        Computation(const std::vector<double>& _input, const Matrix<double>& _W, const std::vector<std::string>& _cellTypes);
        void augmentW(std::vector<std::tuple<std::string,std::string,double>>&);
        std::vector<double> computePerturbation();
        std::vector<double> computeAugmentedPerturbation(); //taking into account virtual nodes in the augmented metapathway

        
};