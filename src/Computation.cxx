

#include "Computation.h"
#include "Matrix.h"
#include "WeightedEdgeGraph.h"
#include "armaUtilities.h"
#include "utilities.h"
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
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

Computation::~Computation(){
    if(metapathway) delete metapathway;
    if(augmentedMetapathway) delete augmentedMetapathway;
}

Computation::Computation(std::string _thisCellType,const std::vector<double>& _input){
    input=_input;
    output=std::vector<double>();
    localCellType = _thisCellType;
    metapathway = new WeightedEdgeGraph();
    augmentedMetapathway = new WeightedEdgeGraph();
    cellTypes = std::vector<std::string>();
}


Computation::Computation(std::string _thisCellType,const std::vector<double>& _input, const Matrix<double>& _W, const std::vector<std::string>& metapathwayNames){
    input=_input;
    output=std::vector<double>();
    localCellType = _thisCellType;
    metapathway = new WeightedEdgeGraph(_W);
    augmentedMetapathway = new WeightedEdgeGraph();
    metapathway->setNodesNames(metapathwayNames); //With default selection of the node names to change(all the nodes in the order established by the matrix rows and columns)
    //augmentedMetapathway = metapathway->addNodes();  // not available since we do not have additional information about the other types
    cellTypes = std::vector<std::string>();


    std::vector<double> normalizationFactors(metapathway->getNumNodes(),0);
    for (int i = 0; i < metapathway->getNumNodes(); i++) {
        for(int j = 0; j < metapathway->getNumNodes();j++){
            normalizationFactors[i] += std::abs(metapathway->getEdgeWeight(i,j)); 
        }
    }
    arma::Mat<double> WtransArma = metapathway->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
    // std::cout<< " WtransArma:\n";
    // print_mat( WtransArma);

    arma::Mat<double> IdentityArma = arma::eye(metapathway->getNumNodes(),metapathway->getNumNodes());
    // std::cout<< "Identity - WtransArma:\n";
    // print_mat(IdentityArma - WtransArma);

    InputArma = Matrix<double>(input).asArmadilloColumnVector();
    pseudoInverseArma = arma::pinv(IdentityArma - WtransArma);
    // std::cout<< "pseudoInverseArma:\n";
    // print_mat(pseudoInverseArma);
    armaInitializedNotAugmented = true;
}

Computation::Computation(std::string _thisCellType,const std::vector<double>& _input, WeightedEdgeGraph* _metapathway, const std::vector<std::string>& metapathwayNames){
    input=_input;
    output=std::vector<double>();
    localCellType = _thisCellType;
    metapathway = _metapathway;
    augmentedMetapathway = new WeightedEdgeGraph();
    metapathway->setNodesNames(metapathwayNames); //With default selection of the node names to change(all the nodes in the order established by the matrix rows and columns)
    //augmentedMetapathway = metapathway->addNodes();  // not available since we do not have additional information about the other types
    cellTypes = std::vector<std::string>();


    std::vector<double> normalizationFactors(metapathway->getNumNodes(),0);
    for (int i = 0; i < metapathway->getNumNodes(); i++) {
        for(int j = 0; j < metapathway->getNumNodes();j++){
            normalizationFactors[i] += std::abs(metapathway->getEdgeWeight(i,j)); 
        }
    }
    arma::Mat<double> WtransArma = metapathway->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
    //TODO normalization by previous weight nodes for the matrix
    
    arma::Mat<double> IdentityArma = arma::eye(metapathway->getNumNodes(),metapathway->getNumNodes());
    InputArma = Matrix<double>(input).asArmadilloColumnVector();
    // std::cout << "[LOG] computing pseudoinverse for metapathway cell : " + localCellType << std::endl;
    // pseudoInverseArma = arma::pinv(IdentityArma - WtransArma);
    // armaInitializedNotAugmented = true;
}

void Computation::augmentMetapathway(const std::vector<std::string>& _celltypes,const std::vector<std::pair<std::string, std::string>>& newEdgesList,const std::vector<double>& newEdgesValue, bool includeSelfVirtual){
    if(augmentedMetapathway) {
        delete augmentedMetapathway;
    }
    try {
        // 
        // std::vector<std::string> tmpcelltypes = _celltypes;
        // auto cellFind = std::find(tmpcelltypes.begin(), tmpcelltypes.end(), localCellType); 
        // if (cellFind != tmpcelltypes.end() && !includeSelfVirtual){
        //     tmpcelltypes.erase(cellFind);  /// PROBLEM!!!
        //     /// this function erases the first element of newEdgesList
        // }
        std::vector<std::string> tmpcelltypes;
        if (!includeSelfVirtual){
            for (uint i = 0; i < _celltypes.size(); i++) {
                if(_celltypes[i] != localCellType) tmpcelltypes.push_back(_celltypes[i]);
            }
            //tmpcelltypes.erase(cellFind);  /// PROBLEM!!!
            /// this function erases the first element of newEdgesList also
        } else {
            tmpcelltypes = _celltypes;
        }
        cellTypes = tmpcelltypes;
        auto virtualNodes = tmpcelltypes;
        for (int i = 0; i < SizeToInt( tmpcelltypes.size()); i++) {
            std::string cellTyp = virtualNodes[i];
            virtualNodes[i] = "v-in:" + cellTyp;
            virtualNodes.push_back("v-out:" + cellTyp);
        }
        augmentedMetapathway = metapathway->addNodesAndCopyNew(virtualNodes);
        for(uint it = 0; it < newEdgesList.size(); it++){
            std::string node1Name = newEdgesList[it].first; 
            std::string node2Name = newEdgesList[it].second;
            double edgeWeight = newEdgesValue[it];

            augmentedMetapathway->addEdge(node1Name,node2Name, edgeWeight);
        }
        std::vector<double> normalizationFactors(augmentedMetapathway->getNumNodes(),0);
        for (int i = 0; i < augmentedMetapathway->getNumNodes(); i++) {
            for(int j = 0; j < augmentedMetapathway->getNumNodes();j++){
                double betaToAdd = std::abs(augmentedMetapathway->getEdgeWeight(i,j));
                normalizationFactors[i] += betaToAdd; 
            }
        }
        arma::Mat<double> WtransAugmentedArma = augmentedMetapathway->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
        //TODO normalization by previous weight nodes for the matrix
        
        arma::Mat<double> IdentityAugmentedArma = arma::eye(augmentedMetapathway->getNumNodes(),augmentedMetapathway->getNumNodes());
        inputAugmented = input;
        for(uint i = 0; i < tmpcelltypes.size()*2; i++){
            inputAugmented.push_back(0.0);
        }
        InputAugmentedArma = arma::Col<double>(inputAugmented);
        std::cout << "[LOG] computing pseudoinverse for augmented metapathway cell : " + localCellType << std::endl;
        pseudoInverseAugmentedArma = arma::pinv(IdentityAugmentedArma - WtransAugmentedArma);
        armaInitializedAugmented = true;


        //get nodeToIndex map as well
        nodeToIndex = augmentedMetapathway->getNodeToIndexMap();
    } catch (...) {
        std::cerr<< "[ERROR] Computation::augmentMetapathway: catch section";
        return;
    }
}

void Computation::augmentMetapathwayNoComputeInverse(const std::vector<std::string>& _celltypes,const std::vector<std::pair<std::string, std::string>>& newEdgesList,const std::vector<double>& newEdgesValue, bool includeSelfVirtual){
    if(augmentedMetapathway) {
        delete augmentedMetapathway;
    }
    try {
        // 
        // std::vector<std::string> tmpcelltypes = _celltypes;
        // auto cellFind = std::find(tmpcelltypes.begin(), tmpcelltypes.end(), localCellType); 
        // if (cellFind != tmpcelltypes.end() && !includeSelfVirtual){
        //     tmpcelltypes.erase(cellFind);  /// PROBLEM!!!
        //     /// this function erases the first element of newEdgesList
        // }
        std::vector<std::string> tmpcelltypes;
        if (!includeSelfVirtual){
            for (uint i = 0; i < _celltypes.size(); i++) {
                if(_celltypes[i] != localCellType) tmpcelltypes.push_back(_celltypes[i]);
            }
            //tmpcelltypes.erase(cellFind);  /// PROBLEM!!!
            /// this function erases the first element of newEdgesList also
        } else {
            tmpcelltypes = _celltypes;
        }
        cellTypes = tmpcelltypes;
        auto virtualNodes = tmpcelltypes;
        for (int i = 0; i < SizeToInt( tmpcelltypes.size()); i++) {
            std::string cellTyp = virtualNodes[i];
            virtualNodes[i] = "v-in:" + cellTyp;
            virtualNodes.push_back("v-out:" + cellTyp);
        }
        augmentedMetapathway = metapathway->addNodesAndCopyNew(virtualNodes);
        for(uint it = 0; it < newEdgesList.size(); it++){
            std::string node1Name = newEdgesList[it].first; 
            std::string node2Name = newEdgesList[it].second;
            double edgeWeight = newEdgesValue[it];

            augmentedMetapathway->addEdge(node1Name,node2Name, edgeWeight);
        }
        std::vector<double> normalizationFactors(augmentedMetapathway->getNumNodes(),0);
        for (int i = 0; i < augmentedMetapathway->getNumNodes(); i++) {
            for(int j = 0; j < augmentedMetapathway->getNumNodes();j++){
                double betaToAdd = std::abs(augmentedMetapathway->getEdgeWeight(i,j));
                normalizationFactors[i] += betaToAdd; 
            }
        }

        // std::cout << "[LOG + TESTING] normalization factors for augmented metapathway cell : " + localCellType << std::endl;
        // printVector(normalizationFactors);
        // std::cout << "[LOG + TESTING] augmented metapathway matrix : " + localCellType << std::endl;
        // augmentedMetapathway->adjMatrix.printMatrix();
        // std::cout << "[LOG + TESTING] augmented metapathway matrix transposed : " + localCellType << std::endl;
        // augmentedMetapathway->adjMatrix.transpose().printMatrix();
        // std::cout << "[LOG + TESTING] augmented metapathway matrix transposed normalized : " + localCellType << std::endl;
        // augmentedMetapathway->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).printMatrix();
        arma::Mat<double> WtransAugmentedArma = augmentedMetapathway->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
        //TODO normalization by previous weight nodes for the matrix
        
        arma::Mat<double> IdentityAugmentedArma = arma::eye(augmentedMetapathway->getNumNodes(),augmentedMetapathway->getNumNodes());
        inputAugmented = input;
        for(uint i = 0; i < tmpcelltypes.size()*2; i++){
            inputAugmented.push_back(0.0);
        }
        InputAugmentedArma = arma::Col<double>(inputAugmented);
        // std::cout << "[LOG] computing pseudoinverse for augmented metapathway cell : " + localCellType << std::endl;
        // pseudoInverseAugmentedArma = arma::pinv(IdentityAugmentedArma - WtransAugmentedArma);
        // armaInitializedAugmented = true;

        //get nodeToIndex map as well
        nodeToIndex = augmentedMetapathway->getNodeToIndexMap();
    } catch (...) {
        std::cerr<< "[ERROR] Computation::augmentMetapathway: catch section";
        return;
    }
}

void Computation::addEdges(const std::vector<std::pair<std::string,std::string>>& newEdgesList, const std::vector<double>& newEdgesValues,bool bothDirections){
    //TODO control over the same length
    int itVal = 0;
    for(auto it = newEdgesList.cbegin(); it!=newEdgesList.cend();it++){
        std::string node1Name = it->first; 
        std::string node2Name = it->second;
        double edgeWeight = newEdgesValues[itVal++];
        if (bothDirections) {
            augmentedMetapathway->addEdge(node2Name,node1Name, edgeWeight);
        }
        augmentedMetapathway->addEdge(node1Name,node2Name, edgeWeight);
    }
    std::vector<double> normalizationFactors(augmentedMetapathway->getNumNodes(),0);
    for (int i = 0; i < augmentedMetapathway->getNumNodes(); i++) {
        for(int j = 0; j < augmentedMetapathway->getNumNodes();j++){
            double betaToAdd = std::abs(augmentedMetapathway->getEdgeWeight(i,j));
            normalizationFactors[i] += betaToAdd; 
        }
    }
    arma::Mat<double> WtransAugmentedArma = augmentedMetapathway->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
    //TODO normalization by previous weight nodes for the matrix
    arma::Mat<double> IdentityAugmentedArma = arma::eye(augmentedMetapathway->getNumNodes(),augmentedMetapathway->getNumNodes());
    std::cout << "[LOG] computing pseudoinverse for augmented metapathway cell : " + localCellType << std::endl;
    pseudoInverseAugmentedArma = arma::pinv(IdentityAugmentedArma - WtransAugmentedArma);
    armaInitializedAugmented = true;
}

void Computation::addEdges(const std::vector<std::tuple<std::string,std::string,double>>& newEdgesList,bool bothDirections){
    //TODO control over the same length
    for(auto it = newEdgesList.cbegin(); it!=newEdgesList.cend();it++){
        std::string node1Name = std::get<0>(*it); 
        std::string node2Name = std::get<1>(*it);
        double edgeWeight = std::get<2>(*it);
        if (bothDirections) {
            augmentedMetapathway->addEdge(node2Name,node1Name, edgeWeight);
        }
        augmentedMetapathway->addEdge(node1Name,node2Name, edgeWeight);
    }
    std::vector<double> normalizationFactors(augmentedMetapathway->getNumNodes(),0);
    for (int i = 0; i < augmentedMetapathway->getNumNodes(); i++) {
        for(int j = 0; j < augmentedMetapathway->getNumNodes();j++){
            double betaToAdd = std::abs(augmentedMetapathway->getEdgeWeight(i,j));
            normalizationFactors[i] += betaToAdd; 
        }
    }
    arma::Mat<double> WtransAugmentedArma = augmentedMetapathway->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
    //TODO normalization by previous weight nodes for the matrix
    arma::Mat<double> IdentityAugmentedArma = arma::eye(augmentedMetapathway->getNumNodes(),augmentedMetapathway->getNumNodes());
    std::cout << "[LOG] computing pseudoinverse for augmented metapathway cell : " + localCellType << std::endl;
    pseudoInverseAugmentedArma = arma::pinv(IdentityAugmentedArma - WtransAugmentedArma);
    armaInitializedAugmented = true;
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

//TODO
std::vector<double> Computation::computeAugmentedPerturbationNorm(){
    arma::Col<double> outputArma =  pseudoInverseAugmentedArma * InputAugmentedArma;
    outputAugmented = armaColumnToVector(outputArma);
    return outputAugmented;
}

std::vector<double> Computation::computeAugmentedPerturbationDissipatedPow2(){
    arma::Col<double> outputArma =  pseudoInverseAugmentedArma * InputAugmentedArma;
    arma::Col<double> dissipationTerm = pow(InputAugmentedArma,2);
    for(uint i = 0;i<outputArma.n_elem;i++){
        outputArma[i] = outputArma[i] - dissipationTerm[i];
    }
    outputAugmented = armaColumnToVector(outputArma);
    return outputAugmented;
}

std::vector<double> Computation::computeAugmentedPerturbationSaturated(const std::vector<double>& saturationsVector){
    if(saturationsVector.size() != 0){
        if (saturationsVector.size() == InputAugmentedArma.n_elem) {
            arma::Col<double> outputArma =  pseudoInverseAugmentedArma * InputAugmentedArma;
            for(uint i = 0;i<outputArma.n_elem;i++){
                outputArma[i] = hyperbolicTangentScaled(outputArma[i], saturationsVector[i]);
            }
            outputAugmented = armaColumnToVector(outputArma);
            return outputAugmented;
        } else{
            throw std::invalid_argument("saturationVector is not of the same size as output vector. abort");
        }
    }
    else {
        arma::Col<double> outputArma =  pseudoInverseAugmentedArma * InputAugmentedArma;
        for(uint i = 0;i<outputArma.n_elem;i++){
            outputArma[i] = hyperbolicTangentScaled(outputArma[i], 1);
        }
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    }
}

std::vector<double> Computation::computeAugmentedPerturbationDissipatedAfterCompute(double timeStep){
    if (dissipationModel) {
        arma::Col<double> outputArma =  pseudoInverseAugmentedArma * InputAugmentedArma;
        arma::Col<double> dissipationTerm = dissipationModel->dissipationTerm(outputArma,timeStep);
        outputArma = outputArma - dissipationTerm;
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    } else {
        throw std::invalid_argument("Computation::computeAugmentedPerturbationDissipatedAfterCompute: dissipationModel is not set");
    }
}

std::vector<double> Computation::computeAugmentedPerturbationDissipatedBeforeCompute(double timeStep){
    if (dissipationModel) {
        arma::Col<double> outputArma =  pseudoInverseAugmentedArma * dissipationModel->dissipate(InputAugmentedArma, timeStep);
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    } else {
        throw std::invalid_argument("Computation::computeAugmentedPerturbationDissipatedBeforeCompute: dissipationModel is not set");
    }
}

std::vector<double> Computation::computeAugmentedPerturbationSaturatedAndDissipatedBeforeCompute(double timeStep, const std::vector<double>& saturationsVector){
    if (saturationsVector.size() ) {
        if (saturationsVector.size() == InputAugmentedArma.n_elem) {
            arma::Col<double> outputArma =  pseudoInverseAugmentedArma * dissipationModel->dissipate(InputAugmentedArma, timeStep);
            for(uint i = 0;i<outputArma.n_elem;i++){
                outputArma[i] = hyperbolicTangentScaled(outputArma[i], saturationsVector[i]);
            }
            outputAugmented = armaColumnToVector(outputArma);
            return outputAugmented;
        } else{
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationSaturatedAndDissipatedBeforeCompute: saturationVector is not of the same size as output vector. abort");
        }
    }
    else {
        arma::Col<double> outputArma =  pseudoInverseAugmentedArma * dissipationModel->dissipate(InputAugmentedArma, timeStep);
        for(uint i = 0;i<outputArma.n_elem;i++){
            outputArma[i] = hyperbolicTangentScaled(outputArma[i], 1);
        }
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    }
}

//TODO
std::vector<double> Computation::computeAugmentedPerturbationEnhanced1(double timeStep, const std::vector<double>& saturationsVector){
    if (saturationsVector.size() ) {
        if (saturationsVector.size() == InputAugmentedArma.n_elem) {
            arma::Col<double> outputArma =  pseudoInverseAugmentedArma * dissipationModel->dissipate(InputAugmentedArma, timeStep);
            for(uint i = 0;i<outputArma.n_elem;i++){
                outputArma[i] = hyperbolicTangentScaled(outputArma[i], saturationsVector[i]);
            }
            outputAugmented = armaColumnToVector(outputArma);
            return outputAugmented;
        } else{
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced1: saturationVector is not of the same size as output vector. abort");
        }
    }
    else {
        arma::Col<double> outputArma =  pseudoInverseAugmentedArma * dissipationModel->dissipate(InputAugmentedArma, timeStep);
        for(uint i = 0;i<outputArma.n_elem;i++){
            outputArma[i] = hyperbolicTangentScaled(outputArma[i], 1);
        }
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    }
}

//TODO
std::vector<double> Computation::computeAugmentedPerturbationEnhanced2(double timeStep, const std::vector<double>& saturationsVector,const std::vector<double>& qVector){
    std::vector<double> saturationVectorVar = saturationsVector;
    std::vector<double> qVectorVar = qVector;
    if(saturationVectorVar.size() == 0){
        saturationVectorVar = std::vector<double>(InputAugmentedArma.n_elem,1);
    }
    if(qVectorVar.size() == 0){
        qVectorVar = std::vector<double>(InputAugmentedArma.n_elem,1);
    }
    if(saturationVectorVar.size() != InputAugmentedArma.n_elem || qVectorVar.size() != InputAugmentedArma.n_elem){
        throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced2: saturationVector or qVector is not of the same size as output vector. abort");
    }
    //dissipation
    arma::Col<double> dissipatedPerturbationArma = dissipationModel->dissipate(InputAugmentedArma, timeStep);
    //saturation
    for(uint i = 0;i<dissipatedPerturbationArma.n_elem;i++){
        dissipatedPerturbationArma[i] = hyperbolicTangentScaled(dissipatedPerturbationArma[i], saturationVectorVar[i]);
    }
    //conservation
    if(augmentedMetapathway == nullptr){
        throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced2: augmentedMetapathway is not set. abort");
    }
    arma::Col<double> outputArma =  dissipatedPerturbationArma - conservationModel->conservationTerm(dissipatedPerturbationArma, augmentedMetapathway->adjMatrix.asArmadilloMatrix(), timeStep);
    outputAugmented = armaColumnToVector(outputArma);
    return outputAugmented;
}


double Computation::getVirtualInputForCell(std::string celltype)const{
    //int index = augmentedMetapathway->getIndexFromName("v-in:" + celltype);
    int index = nodeToIndex.at("v-in:" + celltype);
    if(index > 0) return outputAugmented[index];
    else return 0;

}
double Computation::getVirtualOutputForCell(std::string celltype)const{
    //int index = augmentedMetapathway->getIndexFromName("v-out:" + celltype);
    int index = nodeToIndex.at("v-out:" + celltype);
    if(index > 0) return outputAugmented[index];
    else return 0;
}

void Computation::setInputVinForCell(std::string celltype, double value){
    //TODO if the augmented metapathway is deleted, switch to a direct map saved before deleting the metapathway
    //int index = augmentedMetapathway->getIndexFromName("v-in:" + celltype);
    int index = nodeToIndex.at("v-in:" + celltype);
    if(index > 0) {
        inputAugmented[index]=value;
        InputAugmentedArma[index]=value;
    }
    else throw std::invalid_argument("Computation::setInputVinForCell: invalid set for virtual input: celltype:" + celltype + "does not exist");

}
void Computation::setInputVoutForCell(std::string celltype, double value){
    //int index = augmentedMetapathway->getIndexFromName("v-out:" + celltype);
    int index = nodeToIndex.at("v-out:" + celltype);
    if(index > 0) {
        inputAugmented[index]=value;
        InputAugmentedArma[index]=value;
    }
    else throw std::invalid_argument("Computation::setInputVinForCell: invalid set for virtual input: celltype:" + celltype + "does not exist");
}

void Computation::setDissipationModel(DissipationModel *dissipationModel){
    this->dissipationModel = dissipationModel;
}

void Computation::setConservationModel(ConservationModel *conservationModel){
    this->conservationModel = conservationModel;
}

void Computation::updateInput(const std::vector<double>& newInp, bool augmented){
    if (!augmented) {
        
        if (newInp.size() == 0) {
            InputArma = arma::Col<double>(output);
            input = output;    
        }
        else {
            if(newInp.size() == input.size()){
                InputArma = arma::Col<double>(newInp);
                input = newInp;
            }
        }
    } else {
        if (newInp.size() == 0) {
            InputAugmentedArma = arma::Col<double>(outputAugmented);
            inputAugmented = outputAugmented;    
        }
        else {
            if(newInp.size() == inputAugmented.size()){
                InputAugmentedArma = arma::Col<double>(newInp);
                inputAugmented = newInp;
            }
        }
    }
}


Computation& Computation::operator=( const Computation& rhs){
    metapathway =  rhs.getMetapathway()->copyNew();
    augmentedMetapathway =  rhs.getAugmentedMetapathway()->copyNew();
    input = rhs.getInput();
    output = rhs.getOutput();
    inputAugmented = rhs.getInputAugmented();
    outputAugmented = rhs.getOutputAugmented();
    cellTypes = rhs.getCellTypes();
    localCellType = rhs.getLocalCellType();
    armaInitializedNotAugmented = rhs.isInitializedArmaNotAugmented();
    armaInitializedAugmented = rhs.isInitializedArmaAugmented();
    InputArma = rhs.getInputArma();
    pseudoInverseArma = rhs.getPseudoInverseArma();
    InputAugmentedArma = rhs.getInputAugmentedArma();
    pseudoInverseAugmentedArma = rhs.getPseudoInverseAugmentedArma();
    return *this;
}

Computation Computation::copy()const{
    return *this;
}

void Computation::assign(const Computation & rhs){
    //delete old data in dynamic memory
    if(metapathway){
        delete metapathway;
        metapathway = nullptr;
    }if(augmentedMetapathway){
        delete augmentedMetapathway;
        augmentedMetapathway = nullptr;
    }
    //copy and allocate new structures
    metapathway =  rhs.getMetapathway()->copyNew();
    augmentedMetapathway =  rhs.getAugmentedMetapathway()->copyNew();
    input = rhs.getInput();
    output = rhs.getOutput();
    inputAugmented = rhs.getInputAugmented();
    outputAugmented = rhs.getOutputAugmented();
    cellTypes = rhs.getCellTypes();
    localCellType = rhs.getLocalCellType();
    armaInitializedNotAugmented = rhs.isInitializedArmaNotAugmented();
    armaInitializedAugmented = rhs.isInitializedArmaAugmented();
    InputArma = rhs.getInputArma();
    pseudoInverseArma = rhs.getPseudoInverseArma();
    InputAugmentedArma = rhs.getInputAugmentedArma();
    pseudoInverseAugmentedArma = rhs.getPseudoInverseAugmentedArma();
}

// optimization
void Computation::freeAugmentedGraphs(){
    nodeToIndex = augmentedMetapathway->getNodeToIndexMap();
    delete augmentedMetapathway;
}