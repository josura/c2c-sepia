

#include "Computation.h"
#include "Matrix.h"
#include "WeightedEdgeGraph.h"
#include "armaUtilities.h"
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


    WtransArma = metapathway->adjMatrix.transpose().asArmadilloMatrix();
    IdentityArma = Matrix<double>::createIdentity(metapathway->getNumNodes()).asArmadilloMatrix();
    InputArma = Matrix<double>(input).asArmadilloColumnVector();
    pseudoInverseArma = arma::pinv(IdentityArma - WtransArma);
    armaInitializedNotAugmented = true;
}

void Computation::augmentMetapathway(const std::vector<std::string>& _celltypes,const std::vector<std::tuple<std::string, std::string, double>>& newEdgesList, bool includeSelfVirtual){
    delete augmentedMetapathway;
    try {
        //TODO controls over nodes and edges added? It should be redundand though since the methods of WeightedEdgeGraph are already controlled
        auto cellFind = std::find(_celltypes.begin(), _celltypes.end(), localCellType); 
        std::vector<std::string> tmpcelltypes = _celltypes;
        if (cellFind != _celltypes.end() && !includeSelfVirtual){
            tmpcelltypes.erase(cellFind);
        }
        augmentedMetapathway = metapathway->addNodes(tmpcelltypes);  //these are only one set of nodes
        //TODO differentiate between virtual inputs and virtual outputs
        for(auto it = newEdgesList.cbegin(); it!=newEdgesList.cend();it++){
            std::string node1Name = std::get<0>(*it); 
            std::string node2Name = std::get<1>(*it);
            double edgeWeight = std::get<2>(*it);

            augmentedMetapathway->addEdge(node1Name,node2Name, edgeWeight);
        }
        WtransAugmentedArma = augmentedMetapathway->adjMatrix.transpose().asArmadilloMatrix();
        IdentityAugmentedArma = Matrix<double>::createIdentity(augmentedMetapathway->getNumNodes()).asArmadilloMatrix();
        //TODO augment input with virtual inputs and virtual outputs
        InputAugmentedArma = Matrix<double>(input).asArmadilloColumnVector();
        pseudoInverseAugmentedArma = arma::pinv(IdentityArma - WtransArma);
        armaInitializedAugmented = true;
    } catch (...) {
        std::cerr<< "[ERROR] Computation::augmentMetapathway: catch section";
        return;
    }
}


std::vector<double> Computation::computePerturbation(){
    arma::Col<double> outputArma =  pseudoInverseArma * InputArma;
    return armaColumnToVector(outputArma);
}
std::vector<double> Computation::computeAugmentedPerturbation(){
    arma::Col<double> outputArma =  pseudoInverseAugmentedArma * InputAugmentedArma;
    return armaColumnToVector(outputArma);
}