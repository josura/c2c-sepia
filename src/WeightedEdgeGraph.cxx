#include "WeightedEdgeGraph.h"
#include "utilities.h"
#include <ostream>
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
    this->edgesVector = new std::vector<std::tuple<int, int, float>>();
    
}

WeightedEdgeGraph::WeightedEdgeGraph(int numNodes){
    this->numberOfNodes = numNodes;
    this->nodeWeights = new double[numNodes];
    this->adjList = new std::unordered_set<int>[numNodes];
    this->adjVector = new std::vector<int>[numNodes];
    this->edgesVector = new std::vector<std::tuple<int, int, float>>();
    //edgesVector = new std::vector<std::pair<int, int>>;

    for (int i = 0; i < numNodes; i++) {
        nodeWeights[i]=0;
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

WeightedEdgeGraph* WeightedEdgeGraph::addEdge(int node1, int node2, float weight){
    if(node1 >= numberOfNodes || node2 >= numberOfNodes){
        std::cerr << "add edge failed for edges " << node1 << " and " << node2 << std::endl;
        if(node1 >= numberOfNodes){
            std::cerr << "[ERROR] node1(number "<< node1 << ") is not in the graph that has " << numberOfNodes << " nodes"<<std::endl;
        }
        else{
            std::cerr << "[ERROR] node2(number "<< node2 << ") is not in the graph that has " << numberOfNodes << " nodes"<<std::endl;
        }
    } else if (adjNodes(node1, node2)) {
        //edge already added
    } else {
        numberOfEdges++;
        edgesVector->push_back(std::tuple<int, int, float>(node1,node2, weight));
        adjList[node1].insert(node2);
        adjList[node2].insert(node1);
        adjVector[node1].push_back(node2);
        adjVector[node2].push_back(node1);
    }

    return this;
}


std::tuple<int, int,float>* WeightedEdgeGraph::makeEdgesArray(){
    //emptying old memory
    delete[] edgesArray;
    // new array
    edgesArray = new std::tuple<int, int, float>[numberOfEdges];
    for (int i =0 ; i<numberOfEdges; i++) {
        edgesArray[i] = edgesVector->at(i);
    }
    return edgesArray;
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
        std::cerr << "trying to get an adjacent list of an unknown node: " << node << ">=" << numberOfNodes << std::endl;
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

std::vector<std::tuple<int, int,float>>* WeightedEdgeGraph::getEdgesVector()const{
    return edgesVector;
}

std::tuple<int, int,float>* WeightedEdgeGraph::getEdgesArray()const{
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
    delete edgesVector;

    //creating new data
    this->nodeWeights = new double[g2.numberOfNodes];
    this->adjList = new std::unordered_set<int>[g2.numberOfNodes];
    this->adjVector = new std::vector<int>[g2.numberOfNodes];
    this->edgesVector = new std::vector<std::tuple<int, int, float>>();

    for(auto it = g2.edgesVector->cbegin(); it!=g2.edgesVector->cend();it++){
        this->addEdge(std::get<0>(*it),std::get<1>(*it), std::get<2>(*it));
    }
    makeEdgesArray();
    return *this;
}