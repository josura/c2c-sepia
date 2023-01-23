

#include "Computation.h"
#include "Matrix.h"
#include "WeightedEdgeGraph.h"
#include <iostream>
#include <map>
#include <tuple>
#include <vector>


Computation::Computation(){
    input=std::vector<double>();
    output=std::vector<double>();
    localCellType = "";
    metapathway = new WeightedEdgeGraph();
    augmentedMetapathway = new WeightedEdgeGraph();
    cellTypes = std::vector<std::string>();
}

Computation::Computation(std::string _thisCellType,const std::vector<double>& _input){
    input=_input;
    output=std::vector<double>();
    metapathway = new WeightedEdgeGraph();
    augmentedMetapathway = new WeightedEdgeGraph();
    cellTypes = std::vector<std::string>();
}


Computation::Computation(std::string _thisCellType,const std::vector<double>& _input, const Matrix<double>& _W, const std::vector<std::string>& metapathwayNames){
    input=_input;
    output=std::vector<double>();
    metapathway = new WeightedEdgeGraph(_W);
    metapathway->setNodesNames(metapathwayNames); //With default selection of the node names to change(all the nodes in the order established by the matrix rows and columns)
    //augmentedMetapathway = metapathway->addNodes();  // not available since we do not have additional information about the other types
    cellTypes = std::vector<std::string>();
}

void Computation::augmentMetapathway(const std::vector<std::string>& _celltypes,const std::vector<std::tuple<std::string, std::string, double>>& newEdgesList){
    delete augmentedMetapathway;
    try {
        //TODO controls over nodes and edges added? It should be redundand though since the methods of WeightedEdgeGraph are already controlled
        augmentedMetapathway = metapathway->addNodes(_celltypes);
        for(auto it = newEdgesList.cbegin(); it!=newEdgesList.cend();it++){
            std::string node1Name = std::get<0>(*it); 
            std::string node2Name = std::get<1>(*it);
            double edgeWeight = std::get<2>(*it);

            augmentedMetapathway->addEdge(node1Name,node2Name, edgeWeight);
        }
    } catch (...) {
        std::cerr<< "[ERROR] Computation::augmentMetapathway: catch section";
    }
}