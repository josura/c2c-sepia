

#include "Computation.h"
#include "Matrix.h"
#include "WeightedEdgeGraph.h"
#include <map>
#include <tuple>
#include <vector>


Computation::Computation(){
    input=std::vector<double>();
    output=std::vector<double>();
    metapathway = new WeightedEdgeGraph();
    augmentedMetapathway = new WeightedEdgeGraph();
    geneMapToNode = std::map<std::string,int>();
    cellTypes = std::vector<std::string>();
}

Computation::Computation(const std::vector<double>& _input){
    input=_input;
    output=std::vector<double>();
    metapathway = new WeightedEdgeGraph();
    augmentedMetapathway = new WeightedEdgeGraph();
    geneMapToNode = std::map<std::string,int>();
    cellTypes = std::vector<std::string>();
}

Computation::Computation(const std::vector<double>& _input, const Matrix<double>& _W){
    input=_input;
    output=std::vector<double>();
    metapathway = new WeightedEdgeGraph(_W);
    //augmentedMetapathway = new WeightedEdgeGraph(_W);
    geneMapToNode = std::map<std::string,int>();
    cellTypes = std::vector<std::string>();
}

Computation::Computation(const std::vector<double>& _input,const Matrix<double>& _W,const std::vector<std::string>& _cellTypes){
    input=_input;
    output=std::vector<double>();
    metapathway = new WeightedEdgeGraph(_W);
    //augmentedMetapathway = new WeightedEdgeGraph(_W);
    geneMapToNode = std::map<std::string,int>();
    cellTypes = _cellTypes;
}

void Computation::augmentMetapathway(std::vector<std::string>& _celltypes,std::vector<std::tuple<std::string, std::string, double>>& newEdgesList){
    delete augmentedMetapathway;
    augmentedMetapathway = new WeightedEdgeGraph();
    for(auto it = newEdgesList.cbegin(); it!=newEdgesList.cend();it++){
        std::string node1Name = std::get<0>(*it); 
        std::string node2Name = std::get<1>(*it);
        double edgeWeight = std::get<2>(*it);
        geneMapToNode.insert({node1Name,geneMapToNode.size()});
        geneMapToNode.insert({node2Name,geneMapToNode.size()});

        augmentedMetapathway->addEdge(geneMapToNode[node1Name],geneMapToNode[node2Name], edgeWeight);
    }
    std::vector<std::tuple<std::string, std::string, double>> augmentedEdgesVector;
}