

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
        auto cellFind = std::find(_celltypes.begin(), _celltypes.end(), localCellType); 
        std::vector<std::string> tmpcelltypes = _celltypes;
        if (cellFind != _celltypes.end() && !includeSelfVirtual){
            tmpcelltypes.erase(cellFind);
        }
        auto virtualNodes = tmpcelltypes;
        for (int i = 0; i < tmpcelltypes.size(); i++) {
            virtualNodes[i] = "v-in:" + virtualNodes[i];
            virtualNodes.push_back("v-out:" + virtualNodes[i]);
        }
        augmentedMetapathway = metapathway->addNodes(virtualNodes);
        for(auto it = newEdgesList.cbegin(); it!=newEdgesList.cend();it++){
            std::string node1Name = std::get<0>(*it); 
            std::string node2Name = std::get<1>(*it);
            double edgeWeight = std::get<2>(*it);

            augmentedMetapathway->addEdge(node1Name,node2Name, edgeWeight);
        }
        WtransAugmentedArma = augmentedMetapathway->adjMatrix.transpose().asArmadilloMatrix();
        IdentityAugmentedArma = Matrix<double>::createIdentity(augmentedMetapathway->getNumNodes()).asArmadilloMatrix();
        //TODO augment input(inputAugmented) with virtual inputs and virtual outputs
        inputAugmented = input;
        std::vector<double> zerosVirtualNodes = std::vector<double>(tmpcelltypes.size()*2,0);
        inputAugmented.insert(input.end(),zerosVirtualNodes.begin(),zerosVirtualNodes.end());
        InputAugmentedArma = Matrix<double>(inputAugmented).asArmadilloColumnVector();
        pseudoInverseAugmentedArma = arma::pinv(IdentityArma - WtransArma);
        armaInitializedAugmented = true;
    } catch (...) {
        std::cerr<< "[ERROR] Computation::augmentMetapathway: catch section";
        return;
    }
}


std::vector<double> Computation::computePerturbation(){
    arma::Col<double> outputArma =  pseudoInverseArma * InputArma;
    output = armaColumnToVector(outputArma);
    return output;
}
std::vector<double> Computation::computeAugmentedPerturbation(){
    arma::Col<double> outputArma =  pseudoInverseAugmentedArma * InputAugmentedArma;
    outputAugmented = armaColumnToVector(outputArma);
    return outputAugmented;
}

void Computation::updateInput(const std::vector<double>& newInp, bool augmented){
    if (!augmented) {
        if (newInp.size() == 0) {
            InputArma = Matrix<double>(output).asArmadilloColumnVector();
            input = output;    
        }
        else {
            InputArma = Matrix<double>(newInp).asArmadilloColumnVector();
            input = newInp;
        }
    } else {
        if (newInp.size() == 0) {
            InputAugmentedArma = Matrix<double>(outputAugmented).asArmadilloColumnVector();
            inputAugmented = outputAugmented;    
        }
        else {
            InputAugmentedArma = Matrix<double>(newInp).asArmadilloColumnVector();
            inputAugmented = newInp;
        }
    }
}