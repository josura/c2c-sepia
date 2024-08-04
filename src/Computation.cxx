

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
    graph = new WeightedEdgeGraph();
    augmentedGraph = new WeightedEdgeGraph();
    cellTypes = std::vector<std::string>();
    saturationFunction = [](double value,double saturation)-> double{
        if(value > saturation)return saturation;
        else if (value < -saturation) return -saturation;
        else return value;
    };
}

Computation::~Computation(){
    if(graph) {delete graph; graph=nullptr;}
    if(augmentedGraph) {delete augmentedGraph; augmentedGraph=nullptr;}
}

Computation::Computation(std::string _thisCellType,const std::vector<double>& _input){
    input=_input;
    output=std::vector<double>();
    localCellType = _thisCellType;
    graph = new WeightedEdgeGraph();
    augmentedGraph = new WeightedEdgeGraph();
    cellTypes = std::vector<std::string>();
    saturationFunction = [](double value,double saturation)-> double{
        if(value > saturation)return saturation;
        else if (value < -saturation) return -saturation;
        else return value;
    };
}


Computation::Computation(std::string _thisCellType,const std::vector<double>& _input, const Matrix<double>& _W, const std::vector<std::string>& graphNames){
    input=_input;
    output=std::vector<double>();
    localCellType = _thisCellType;
    graph = new WeightedEdgeGraph(_W);
    augmentedGraph = new WeightedEdgeGraph();
    graph->setNodesNames(graphNames); //With default selection of the node names to change(all the nodes in the order established by the matrix rows and columns)
    cellTypes = std::vector<std::string>();


    std::vector<double> normalizationFactors(graph->getNumNodes(),0);
    for (int i = 0; i < graph->getNumNodes(); i++) {
        for(int j = 0; j < graph->getNumNodes();j++){
            normalizationFactors[i] += std::abs(graph->getEdgeWeight(i,j)); 
        }
    }
    arma::Mat<double> WtransArma = graph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
    
    arma::Mat<double> IdentityArma = arma::eye(graph->getNumNodes(),graph->getNumNodes());
    
    InputArma = Matrix<double>(input).asArmadilloColumnVector();
    pseudoInverseArma = arma::pinv(IdentityArma - WtransArma);
    armaInitializedNotAugmented = true;
    saturationFunction = [](double value,double saturation)-> double{
        if(value > saturation)return saturation;
        else if (value < -saturation) return -saturation;
        else return value;
    };
}

Computation::Computation(std::string _thisCellType,const std::vector<double>& _input, WeightedEdgeGraph* _graph, const std::vector<std::string>& graphNames){
    input=_input;
    output=std::vector<double>();
    localCellType = _thisCellType;
    graph = _graph;
    augmentedGraph = new WeightedEdgeGraph();
    graph->setNodesNames(graphNames); //With default selection of the node names to change(all the nodes in the order established by the matrix rows and columns)
    cellTypes = std::vector<std::string>();


    std::vector<double> normalizationFactors(graph->getNumNodes(),0);
    for (int i = 0; i < graph->getNumNodes(); i++) {
        for(int j = 0; j < graph->getNumNodes();j++){
            normalizationFactors[i] += std::abs(graph->getEdgeWeight(i,j)); 
        }
    }
    arma::Mat<double> WtransArma = graph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
    
    arma::Mat<double> IdentityArma = arma::eye(graph->getNumNodes(),graph->getNumNodes());
    InputArma = Matrix<double>(input).asArmadilloColumnVector();
    // std::cout << "[LOG] computing pseudoinverse for graph cell : " + localCellType << std::endl;
    // pseudoInverseArma = arma::pinv(IdentityArma - WtransArma);
    // armaInitializedNotAugmented = true;
    saturationFunction = [](double value,double saturation)-> double{
        if(value > saturation)return saturation;
        else if (value < -saturation) return -saturation;
        else return value;
    };
}

void Computation::augmentGraph(const std::vector<std::string>& _types,const std::vector<std::pair<std::string, std::string>>& newEdgesList,const std::vector<double>& newEdgesValue, bool includeSelfVirtual){
    if(augmentedGraph) {
        delete augmentedGraph;
    }
    try {
        std::vector<std::string> tmptypes;
        if (!includeSelfVirtual){
            for (uint i = 0; i < _types.size(); i++) {
                if(_types[i] != localCellType) tmptypes.push_back(_types[i]);
            }
        } else {
            tmptypes = _types;
        }
        cellTypes = tmptypes;
        auto virtualNodes = tmptypes;
        for (int i = 0; i < SizeToInt( tmptypes.size()); i++) {
            std::string cellTyp = virtualNodes[i];
            virtualNodes[i] = "v-in:" + cellTyp;
            virtualNodes.push_back("v-out:" + cellTyp);
        }
        augmentedGraph = graph->addNodesAndCopyNew(virtualNodes);
        for(uint it = 0; it < newEdgesList.size(); it++){
            std::string node1Name = newEdgesList[it].first; 
            std::string node2Name = newEdgesList[it].second;
            double edgeWeight = newEdgesValue[it];

            augmentedGraph->addEdge(node1Name,node2Name, edgeWeight);
        }
        std::vector<double> normalizationFactors(augmentedGraph->getNumNodes(),0);
        for (int i = 0; i < augmentedGraph->getNumNodes(); i++) {
            for(int j = 0; j < augmentedGraph->getNumNodes();j++){
                double betaToAdd = std::abs(augmentedGraph->getEdgeWeight(i,j));
                normalizationFactors[i] += betaToAdd; 
            }
        }
        arma::Mat<double> WtransAugmentedArma = augmentedGraph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
        //TODO normalization by previous weight nodes for the matrix
        
        arma::Mat<double> IdentityAugmentedArma = arma::eye(augmentedGraph->getNumNodes(),augmentedGraph->getNumNodes());
        inputAugmented = input;
        for(uint i = 0; i < tmptypes.size()*2; i++){
            inputAugmented.push_back(0.0);
        }
        InputAugmentedArma = arma::Col<double>(inputAugmented);
        std::cout << "[LOG] computing pseudoinverse for augmented graph cell : " + localCellType << std::endl;
        pseudoInverseAugmentedArma = arma::pinv(IdentityAugmentedArma - WtransAugmentedArma);
        armaInitializedAugmented = true;


        //get nodeToIndex map as well
        nodeToIndex = augmentedGraph->getNodeToIndexMap();
    } catch (...) {
        std::cerr<< "[ERROR] Computation::augmentGraph: catch section";
        return;
    }
}

void Computation::augmentGraphNoComputeInverse(const std::vector<std::string>& _types,const std::vector<std::pair<std::string, std::string>>& newEdgesList,const std::vector<double>& newEdgesValue, bool includeSelfVirtual){
    if(augmentedGraph) {
        delete augmentedGraph;
    }
    try {
        std::vector<std::string> tmptypes;
        if (!includeSelfVirtual){
            for (uint i = 0; i < _types.size(); i++) {
                if(_types[i] != localCellType) tmptypes.push_back(_types[i]);
            }
        } else {
            tmptypes = _types;
        }
        cellTypes = tmptypes;
        auto virtualNodes = tmptypes;
        for (int i = 0; i < SizeToInt( tmptypes.size()); i++) {
            std::string cellTyp = virtualNodes[i];
            virtualNodes[i] = "v-in:" + cellTyp;
            virtualNodes.push_back("v-out:" + cellTyp);
        }
        augmentedGraph = graph->addNodesAndCopyNew(virtualNodes);
        for(uint it = 0; it < newEdgesList.size(); it++){
            std::string node1Name = newEdgesList[it].first; 
            std::string node2Name = newEdgesList[it].second;
            double edgeWeight = newEdgesValue[it];

            augmentedGraph->addEdge(node1Name,node2Name, edgeWeight);
        }
        std::vector<double> normalizationFactors(augmentedGraph->getNumNodes(),0);
        for (int i = 0; i < augmentedGraph->getNumNodes(); i++) {
            for(int j = 0; j < augmentedGraph->getNumNodes();j++){
                double betaToAdd = std::abs(augmentedGraph->getEdgeWeight(i,j));
                normalizationFactors[i] += betaToAdd; 
            }
        }

        // std::cout << "[LOG + TESTING] normalization factors for augmented graph cell : " + localCellType << std::endl;
        // printVector(normalizationFactors);
        // std::cout << "[LOG + TESTING] augmented graph matrix : " + localCellType << std::endl;
        // augmentedGraph->adjMatrix.printMatrix();
        // std::cout << "[LOG + TESTING] augmented graph matrix transposed : " + localCellType << std::endl;
        // augmentedGraph->adjMatrix.transpose().printMatrix();
        // std::cout << "[LOG + TESTING] augmented graph matrix transposed normalized : " + localCellType << std::endl;
        // augmentedGraph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).printMatrix();
        arma::Mat<double> WtransAugmentedArma = augmentedGraph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
        //TODO normalization by previous weight nodes for the matrix
        
        arma::Mat<double> IdentityAugmentedArma = arma::eye(augmentedGraph->getNumNodes(),augmentedGraph->getNumNodes());
        inputAugmented = input;
        for(uint i = 0; i < tmptypes.size()*2; i++){
            inputAugmented.push_back(0.0);
        }
        InputAugmentedArma = arma::Col<double>(inputAugmented);
        // std::cout << "[LOG] computing pseudoinverse for augmented graph cell : " + localCellType << std::endl;
        // pseudoInverseAugmentedArma = arma::pinv(IdentityAugmentedArma - WtransAugmentedArma);
        // armaInitializedAugmented = true;

        //get nodeToIndex map as well
        nodeToIndex = augmentedGraph->getNodeToIndexMap();
    } catch (...) {
        std::cerr<< "[ERROR] Computation::augmentGraph: catch section";
        return;
    }
}

void Computation::addEdges(const std::vector<std::pair<std::string,std::string>>& newEdgesList, const std::vector<double>& newEdgesValues,bool bothDirections, bool inverseComputation){
    //TODO control over the same length
    int itVal = 0;
    for(auto it = newEdgesList.cbegin(); it!=newEdgesList.cend();it++){
        std::string node1Name = it->first; 
        std::string node2Name = it->second;
        double edgeWeight = newEdgesValues[itVal++];
        if (bothDirections) {
            augmentedGraph->addEdge(node2Name,node1Name, edgeWeight);
        }
        augmentedGraph->addEdge(node1Name,node2Name, edgeWeight);
    }
    std::vector<double> normalizationFactors(augmentedGraph->getNumNodes(),0);
    for (int i = 0; i < augmentedGraph->getNumNodes(); i++) {
        for(int j = 0; j < augmentedGraph->getNumNodes();j++){
            double betaToAdd = std::abs(augmentedGraph->getEdgeWeight(i,j));
            normalizationFactors[i] += betaToAdd; 
        }
    }
    if(inverseComputation){
        arma::Mat<double> WtransAugmentedArma = augmentedGraph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
        //TODO normalization by previous weight nodes for the matrix
        arma::Mat<double> IdentityAugmentedArma = arma::eye(augmentedGraph->getNumNodes(),augmentedGraph->getNumNodes());
        std::cout << "[LOG] computing pseudoinverse for augmented graph cell : " + localCellType << std::endl;
        pseudoInverseAugmentedArma = arma::pinv(IdentityAugmentedArma - WtransAugmentedArma);
        armaInitializedAugmented = true;
    }
    
}

void Computation::addEdges(const std::vector<std::tuple<std::string,std::string,double>>& newEdgesList,bool bothDirections, bool inverseComputation){
    //TODO control over the same length
    for(auto it = newEdgesList.cbegin(); it!=newEdgesList.cend();it++){
        std::string node1Name = std::get<0>(*it); 
        std::string node2Name = std::get<1>(*it);
        double edgeWeight = std::get<2>(*it);
        if (bothDirections) {
            augmentedGraph->addEdge(node2Name,node1Name, edgeWeight);
        }
        augmentedGraph->addEdge(node1Name,node2Name, edgeWeight);
    }
    std::vector<double> normalizationFactors(augmentedGraph->getNumNodes(),0);
    for (int i = 0; i < augmentedGraph->getNumNodes(); i++) {
        for(int j = 0; j < augmentedGraph->getNumNodes();j++){
            double betaToAdd = std::abs(augmentedGraph->getEdgeWeight(i,j));
            normalizationFactors[i] += betaToAdd; 
        }
    }
    if(inverseComputation){
        arma::Mat<double> WtransAugmentedArma = augmentedGraph->adjMatrix.transpose().normalizeByVectorColumn(normalizationFactors).asArmadilloMatrix();
        //TODO normalization by previous weight nodes for the matrix
        arma::Mat<double> IdentityAugmentedArma = arma::eye(augmentedGraph->getNumNodes(),augmentedGraph->getNumNodes());
        std::cout << "[LOG] computing pseudoinverse for augmented graph cell : " + localCellType << std::endl;
        pseudoInverseAugmentedArma = arma::pinv(IdentityAugmentedArma - WtransAugmentedArma);
        armaInitializedAugmented = true;
    }
}


void Computation::addEdgesAndNodes(const std::vector<std::tuple<std::string,std::string,double>>& newEdgesList,bool bothDirections, bool inverseComputation){
    // get the nodes that are not in the graph yet and add them
    std::vector<std::string> nodesToAdd;
    for(auto it = newEdgesList.cbegin(); it != newEdgesList.cend(); it++){
        std::string node1Name = std::get<0>(*it); 
        std::string node2Name = std::get<1>(*it);        

        if(!augmentedGraph->containsNode(node1Name) && !vectorContains(nodesToAdd,node1Name)){
            nodesToAdd.push_back(node1Name);
        }
        if(!augmentedGraph->containsNode(node2Name) && !vectorContains(nodesToAdd,node2Name)){
            nodesToAdd.push_back(node2Name);
        }
    }
    augmentedGraph->addNodes(nodesToAdd);

    //get nodeToIndex map as well
    nodeToIndex = augmentedGraph->getNodeToIndexMap();
    // update inputAugmented and InputAugmentedArma
    for(uint i = 0; i < nodesToAdd.size(); i++){
        inputAugmented.push_back(0.0);
    }
    InputAugmentedArma = arma::Col<double>(inputAugmented);
    // add the edges
    this->addEdges(newEdgesList,bothDirections,inverseComputation);

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

//TODO test if it works in the correct way
std::vector<double> Computation::computeAugmentedPerturbationEnhanced1(double timeStep, bool saturation, const std::vector<double>& saturationsVector){
    if(saturation){
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
    } else {
        arma::Col<double> outputArma =  pseudoInverseAugmentedArma * dissipationModel->dissipate(InputAugmentedArma, timeStep);
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    }
}

//TODO test if it works in the correct way
std::vector<double> Computation::computeAugmentedPerturbationEnhanced2(double timeStep, bool saturation, const std::vector<double>& saturationsVector,const std::vector<double>& qVector){
    if (saturation) {
        std::vector<double> saturationVectorVar = saturationsVector;
        std::vector<double> qVectorVar = qVector;
        if(saturationVectorVar.size() == 0){
            saturationVectorVar = std::vector<double>(InputAugmentedArma.n_elem,1);
        }
        if(saturationVectorVar.size() != InputAugmentedArma.n_elem ){
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced2: saturationVector is not of the same size as output vector. abort");
        }
        //dissipation
        arma::Col<double> dissipatedPerturbationArma = dissipationModel->dissipate(InputAugmentedArma, timeStep);
        //conservation
        if(augmentedGraph == nullptr){
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced2: augmentedGraph is not set. abort");
        }
        arma::Col<double> outputArma =  pseudoInverseAugmentedArma * dissipatedPerturbationArma - conservationModel->conservationTerm(dissipatedPerturbationArma, normalize1Rows(augmentedGraph->adjMatrix.asArmadilloMatrix()) , timeStep, qVectorVar);
        //saturation
        for(uint i = 0;i<outputArma.n_elem;i++){
            double saturatedValue = hyperbolicTangentScaled(outputArma[i], saturationVectorVar[i]);
            outputArma[i] = saturatedValue;
        }
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    } else {
        //dissipation
        arma::Col<double> dissipatedPerturbationArma = dissipationModel->dissipate(InputAugmentedArma, timeStep);
        //conservation
        if(augmentedGraph == nullptr){
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced2: augmentedGraph is not set. abort");
        }
        arma::Col<double> conservationVector = conservationModel->conservationTerm(dissipatedPerturbationArma, normalize1Rows(augmentedGraph->adjMatrix.asArmadilloMatrix()), timeStep, qVector);
        arma::Col<double> outputArma =  pseudoInverseAugmentedArma * dissipatedPerturbationArma - conservationVector;
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    }
}

std::vector<double> Computation::computeAugmentedPerturbationEnhanced3(double timeStep, bool saturation, const std::vector<double>& saturationsVector, const std::vector<double>& qVector, std::function<double(double)> propagationScaleFunction){
    if (saturation) {
        std::vector<double> saturationVectorVar = saturationsVector;
        std::vector<double> qVectorVar = qVector;
        if(saturationVectorVar.size() == 0){
            saturationVectorVar = std::vector<double>(InputAugmentedArma.n_elem,1);
        }
        if(saturationVectorVar.size() != InputAugmentedArma.n_elem ){
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced3: saturationVector is not of the same size as output vector. abort");
        }
        //dissipation
        arma::Col<double> dissipatedPerturbationArma = dissipationModel->dissipate(InputAugmentedArma, timeStep);
        //conservation
        if(augmentedGraph == nullptr){
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced3: augmentedGraph is not set. abort");
        }
        arma::Col<double> outputArma = pseudoInverseAugmentedArma * dissipatedPerturbationArma * propagationScaleFunction(timeStep) - conservationModel->conservationTerm(dissipatedPerturbationArma, normalize1Rows(augmentedGraph->adjMatrix.asArmadilloMatrix()) , timeStep, qVectorVar);
        //saturation
        for(uint i = 0;i<outputArma.n_elem;i++){
            double saturatedValue = hyperbolicTangentScaled(outputArma[i], saturationVectorVar[i]);
            outputArma[i] = saturatedValue;
        }
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    } else {
        //dissipation
        arma::Col<double> dissipatedPerturbationArma = dissipationModel->dissipate(InputAugmentedArma, timeStep);
        //conservation
        if(augmentedGraph == nullptr){
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced3: augmentedGraph is not set. abort");
        }
        arma::Col<double> conservationVector = conservationModel->conservationTerm(dissipatedPerturbationArma, normalize1Rows(augmentedGraph->adjMatrix.asArmadilloMatrix()), timeStep, qVector);
        arma::Col<double> outputArma = pseudoInverseAugmentedArma * dissipatedPerturbationArma * propagationScaleFunction(timeStep) - conservationVector;
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    }
}


std::vector<double> Computation::computeAugmentedPerturbationEnhanced4(double timeStep, bool saturation, const std::vector<double>& saturationsVector, const std::vector<double>& qVector){
    if (saturation) {
        std::vector<double> saturationVectorVar = saturationsVector;
        std::vector<double> qVectorVar = qVector;
        if(saturationVectorVar.size() == 0){
            saturationVectorVar = std::vector<double>(InputAugmentedArma.n_elem,1);
        }
        if(saturationVectorVar.size() != InputAugmentedArma.n_elem ){
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced4: saturationVector is not of the same size as output vector. abort");
        }
        //dissipation
        arma::Col<double> dissipatedPerturbationArma;
        try
        {
            dissipatedPerturbationArma = dissipationModel->dissipate(InputAugmentedArma, timeStep);
            //conservation
            if(augmentedGraph == nullptr){
                throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced4: augmentedGraph is not set. abort");
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced4: error during the computation of dissipation");
        }
        arma::Col<double> outputArma = propagationModel->propagate(dissipatedPerturbationArma,timeStep) - conservationModel->conservationTerm(dissipatedPerturbationArma, normalize1Rows(augmentedGraph->adjMatrix.asArmadilloMatrix()) , timeStep, qVectorVar);
        //saturation
        outputAugmented = armaColumnToVector(outputArma);
        for(uint i = 0;i<outputAugmented.size();i++){
            // TODO change to be more general, doesn't make sense to have the hyperbolic tangent if the condition is controlled (already a control of saturation)
            if(std::abs(outputAugmented[i]) > saturationVectorVar[i]){
                //double saturatedValue = hyperbolicTangentScaled(outputAugmented[i], saturationVectorVar[i]); //problems when the ouput value is big enough, gives nan
                double saturatedValue = saturationFunction(outputAugmented[i], saturationVectorVar[i]);
                outputAugmented[i] = saturatedValue;
            }
        }
        return outputAugmented;
    } else {
        //dissipation
        arma::Col<double> dissipatedPerturbationArma = dissipationModel->dissipate(InputAugmentedArma, timeStep);
        //conservation
        if(augmentedGraph == nullptr){
            throw std::invalid_argument("[ERROR] Computation::computeAugmentedPerturbationEnhanced4: augmentedGraph is not set. abort");
        }
        arma::Col<double> conservationVector = conservationModel->conservationTerm(dissipatedPerturbationArma, normalize1Rows(augmentedGraph->adjMatrix.asArmadilloMatrix()), timeStep, qVector);
        arma::Col<double> outputArma = propagationModel->propagate(dissipatedPerturbationArma,timeStep) - conservationVector;
        outputAugmented = armaColumnToVector(outputArma);
        return outputAugmented;
    }
}

        

double Computation::getVirtualInputForType(std::string type, std::string sourceNode)const{
    if(sourceNode != ""){
        type = type + "_" + sourceNode;
    }
    int index = nodeToIndex.at("v-in:" + type);
    if(index > 0) return outputAugmented[index];
    else return 0;

}
double Computation::getVirtualOutputForType(std::string type, std::string targetNode)const{
    if(targetNode != ""){
        type = type + "_" + targetNode;
    }
    int index = nodeToIndex.at("v-out:" + type);
    if(index > 0) return outputAugmented[index];
    else return 0;
}

void Computation::setInputVinForType(std::string type, double value, std::string sourceNode){
    //TODO if the augmented graph is deleted, switch to a direct map saved before deleting the graph
    //int index = augmentedGraph->getIndexFromName("v-in:" + type);
    std::string vinName = "v-in:" + type;
    if(sourceNode != ""){
        vinName = vinName + "_" + sourceNode;
    }
    int index = nodeToIndex.at("v-in:" + type);
    if(index > 0) {
        inputAugmented[index]=value;
        InputAugmentedArma[index]=value;
    }
    else throw std::invalid_argument("Computation::setInputVinForType: invalid set for virtual input: type:" + type + "does not exist");

}
void Computation::setInputVoutForType(std::string type, double value, std::string targetNode){
    //int index = augmentedGraph->getIndexFromName("v-out:" + type);
    std::string voutName = "v-out:" + type;
    if(targetNode != ""){
        voutName = voutName + "_" + targetNode;
    }
    int index = nodeToIndex.at(voutName);
    if(index > 0) {
        inputAugmented[index]=value;
        InputAugmentedArma[index]=value;
    }
    else throw std::invalid_argument("Computation::setInputVinForType: invalid set for virtual input: type:" + type + "does not exist");
}

void Computation::setDissipationModel(DissipationModel *dissipationModel){
    this->dissipationModel = dissipationModel;
}

void Computation::setConservationModel(ConservationModel *conservationModel){
    this->conservationModel = conservationModel;
}

void Computation::setPropagationModel(PropagationModel *propagationModel){
    this->propagationModel = propagationModel;
}


void Computation::setInputAugmented(const std::vector<double>& newInputAugmented){
    if(newInputAugmented.size() == inputAugmented.size()){
        inputAugmented = newInputAugmented;
        InputAugmentedArma = arma::Col<double>(newInputAugmented);
    }
}

void Computation::resetVirtualOutputs(){
    //virtual outputs have the names starting with v-out:
    for(auto it = nodeToIndex.cbegin(); it!=nodeToIndex.cend();it++){
        if(it->first.find("v-out:") != std::string::npos){
            InputAugmentedArma[it->second] = 0;
        }
    } 

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
    graph =  rhs.getGraph()->copyNew();
    augmentedGraph =  rhs.getAugmentedGraph()->copyNew();
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
    if(graph){
        delete graph;
        graph = nullptr;
    }if(augmentedGraph){
        delete augmentedGraph;
        augmentedGraph = nullptr;
    }
    //copy and allocate new structures
    graph =  rhs.getGraph()->copyNew();
    augmentedGraph =  rhs.getAugmentedGraph()->copyNew();
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
    nodeToIndex = augmentedGraph->getNodeToIndexMap();
    delete augmentedGraph;
}