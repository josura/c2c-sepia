#include "WeightedEdgeGraph.h"
#include "Matrix.h"
#include "utilities.h"
//#include "utilities.h"
#include <ostream>
#include <stdexcept>
#include <string>
#include <sys/types.h>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <utility>


WeightedEdgeGraph::WeightedEdgeGraph(){
    this->numberOfNodes = 0;
    this->nodeWeights = nullptr;
    this->adjList = nullptr;
    this->adjVector = nullptr;
    this->edgesArray = nullptr;
    this->adjMatrix = Matrix<double>();
    
}

WeightedEdgeGraph::WeightedEdgeGraph(int numNodes){
    this->numberOfNodes = numNodes;
    this->nodeWeights = new double[numNodes];
    this->adjList = new std::unordered_set<int>[numNodes];
    this->adjVector = new std::vector<int>[numNodes];
    //edgesVector = new std::vector<std::pair<int, int>>;

    for (int i = 0; i < numNodes; i++) {
        nodeWeights[i]=0;
    }
    this->adjMatrix = Matrix<double>(numNodes,numNodes);
    
}


WeightedEdgeGraph::WeightedEdgeGraph(const Matrix<double>& _adjMatrix){
    if (_adjMatrix.getCols()==_adjMatrix.getRows()) {
        int numNodes = _adjMatrix.getCols();
        this->numberOfNodes = numNodes;
        this->nodeWeights = new double[numNodes];
        for (int i = 0; i < numNodes; i++) {
            nodeWeights[i]=0;
        }
        this->adjList = new std::unordered_set<int>[numNodes];
        this->adjVector = new std::vector<int>[numNodes];
        this->adjMatrix = Matrix<double>(numNodes,numNodes);
        for (int i = 0 ; i<numNodes; i++) {
            for (int j = 0; j<numNodes; j++) {
                if(!approximatelyEqual(_adjMatrix.getValue(i, j),0.0,0.0000000001)) 
                    addEdge(i, j, _adjMatrix.getValue(i, j));
            }
        }
    }
    
}

WeightedEdgeGraph::~WeightedEdgeGraph(){
    //for (int i=0; i<numberOfNodes; i++) {
    //    delete this->adjList[i];
    //}
    delete [] this->adjList;
    delete [] this->nodeWeights;
    
}

int WeightedEdgeGraph::degreeOfNode(int node)const{
    return adjList[node].size();
}

WeightedEdgeGraph* WeightedEdgeGraph::addEdge(int node1, int node2, double weight){
    if(node1 >= numberOfNodes || node2 >= numberOfNodes){
        std::cerr << "add edge failed for edges " << std::to_string(node1) << " and " << std::to_string(node2) << std::endl;
        if(node1 >= numberOfNodes){
            std::cerr << "[ERROR] node1(number "<< std::to_string(node1) << ") is not in the graph that has " << std::to_string(numberOfNodes) << " nodes"<<std::endl;
        }
        else{
            std::cerr << "[ERROR] node2(number "<< std::to_string(node2) << ") is not in the graph that has " << std::to_string(numberOfNodes) << " nodes"<<std::endl;
        }
    } else if (adjNodes(node1, node2)) {
        //edge already added
    } else {
        numberOfEdges++;
        edgesVector.push_back(std::tuple<int, int, double>(node1,node2, weight));
        adjList[node1].insert(node2);
        adjList[node2].insert(node1);
        adjVector[node1].push_back(node2);
        adjVector[node2].push_back(node1);
        adjMatrix(node1,node2) = weight;
    }

    return this;
}


WeightedEdgeGraph* WeightedEdgeGraph::addEdge(std::string node1, std::string node2, double weight){

}

WeightedEdgeGraph* WeightedEdgeGraph::addNode(double value){

}
WeightedEdgeGraph* WeightedEdgeGraph::addNode(std::string name, double value){

}


WeightedEdgeGraph* WeightedEdgeGraph::addNodes(std::vector<double> values){

}
WeightedEdgeGraph* WeightedEdgeGraph::addNodes(std::vector<std::string> names, std::vector<double> values){
    
}

WeightedEdgeGraph* WeightedEdgeGraph::setNodeValue(int node, double value){

}
WeightedEdgeGraph* WeightedEdgeGraph::setNodeValue(std::string node, double value){

}


std::tuple<int, int,double>* WeightedEdgeGraph::makeEdgesArray(){
    //emptying old memory
    delete[] edgesArray;
    // new array
    edgesArray = new std::tuple<int, int, double>[numberOfEdges];
    for (int i =0 ; i<numberOfEdges; i++) {
        edgesArray[i] = edgesVector.at(i);
    }
    return edgesArray;
}

//TODO complete
Matrix<double> WeightedEdgeGraph::makeMatrix(){
    //emptying old memory
    // new array
    return this->adjMatrix;
}


// accessory functions

int WeightedEdgeGraph::getNumNodes()const {
    return numberOfNodes;
}

int WeightedEdgeGraph::getNumEdges()const {
    return numberOfEdges;
}

double* WeightedEdgeGraph::getNodeWeights()const{
    return nodeWeights;
}

std::string WeightedEdgeGraph::getNodeWeightsStr()const{
    std::string stringa = "";
    for (int i = 0; i<numberOfNodes; i++) {
        stringa += std::to_string(nodeWeights[i]) + std::string(" ");
    }
    return stringa;
}


std::unordered_set<int>* WeightedEdgeGraph::getAdjList(int node)const{
    if(node>=numberOfNodes){
        std::cerr << "trying to get an adjacent list of an unknown node: " << std::to_string(node) << ">=" << std::to_string(numberOfNodes) << std::endl;
        return NULL;
    }
    return &adjList[node];
}

std::string WeightedEdgeGraph::getAdjListStr(int node)const{
    std::string stringa;
    for(auto it = adjList[node].cbegin(); it != adjList[node].cend();it++){
                stringa += std::to_string(*it) + " ";
            }
    return stringa;
}

std::vector<std::tuple<int, int,double>> WeightedEdgeGraph::getEdgesVector()const{
    return edgesVector;
}

std::tuple<int, int,double>* WeightedEdgeGraph::getEdgesArray()const{
    return edgesArray;
}

bool WeightedEdgeGraph::adjNodes(int node1, int node2){
    return ( (adjList[node1].find(node2) != adjList[node1].end()) || (adjList[node2].find(node1) != adjList[node2].end())) ;
}

WeightedEdgeGraph& WeightedEdgeGraph::operator=(const WeightedEdgeGraph& g2){
    this->numberOfNodes = g2.numberOfNodes;
    //deleting old data
    delete[] nodeWeights;
    delete[] adjList;
    delete[] adjVector;
    delete[] edgesArray;

    //creating new data
    this->nodeWeights = new double[g2.numberOfNodes];
    this->adjList = new std::unordered_set<int>[g2.numberOfNodes];
    this->adjVector = new std::vector<int>[g2.numberOfNodes];
    this->edgesVector.clear();
    this->adjMatrix = Matrix<double>(g2.numberOfNodes,g2.numberOfNodes);

    for(auto it = g2.edgesVector.cbegin(); it!=g2.edgesVector.cend();it++){
        this->addEdge(std::get<0>(*it),std::get<1>(*it), std::get<2>(*it));
    }
    makeEdgesArray();
    return *this;
}

std::ostream& operator<< (std::ostream &out, const WeightedEdgeGraph& data) {
            out << data.getNumNodes() << " " << data.getNumEdges() <<std::endl;
            std::string nodeweights = data.getNodeWeightsStr();
            out << nodeweights << std::endl;
            out << "Adj Lists" << std::endl;
            for(int i = 0; i<data.getNumNodes(); i++){
                out << "node " << i << " :" << data.getAdjListStr(i) << std::endl;
            }
            out << "Edges vector: {";
            std::vector<std::tuple<int,int,double>> tmpVec= data.getEdgesVector();
            for(auto it = tmpVec.cbegin();it!=tmpVec.cend();it++){
                out << "(" << get<0>(*it)<< ","<< get<1>(*it)<< ","<< get<2>(*it)<< "," << ")," ;
            }
            out << "}"<< std::endl;
            return out;
        }


// WeightedEdgeGraph& WeightedEdgeGraph::setAdjMatrix(Matrix<double>& mat){
//     if (mat.getCols()==mat.getRows()) {
//         this->numberOfNodes = mat.getCols();
//         //deleting old data
//         delete[] nodeWeights;
//         delete[] adjList;
//         delete[] adjVector;
//         delete[] edgesArray;

//         //creating new data
//         this->nodeWeights = new double[numberOfNodes];
//         this->adjList = new std::unordered_set<int>[this->numberOfNodes];
//         this->adjVector = new std::vector<int>[this->numberOfNodes];
//         this->edgesVector.clear();
//         for (int i = 0 ; i<this->numberOfNodes; i++) {
//             for (int j = 0; j<this->numberOfNodes; j++) {
//                 addEdge(i, j, mat.getValue(i, j));
//             }
//         }
//         return *this;
//     } else throw std::invalid_argument("adjMatrix is not square(does not represent a graph)");
// }