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
    numberOfNodes = 0;
    this->nodeWeights = nullptr;
    adjList = nullptr;
    adjVector = nullptr;
    
}

WeightedEdgeGraph::WeightedEdgeGraph(int numNodes){
    numberOfNodes = numNodes;
    this->nodeWeights = new double[numNodes];
    adjList = new std::unordered_set<int>[numNodes];
    adjVector = new std::vector<int>[numNodes];
    //edgesVector = new std::vector<std::pair<int, int>>;
    
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
    } else if (adjNodes(node1, node2)) {
        //edge already added
    } else {
        numberOfEdges++;
        edgesVector.push_back(std::tuple<int, int, float>(node1,node2, weight));
        adjList[node1].insert(node2);
        adjList[node2].insert(node1);
        adjVector[node1].push_back(node2);
        adjVector[node2].push_back(node1);
    }

    return this;
}


std::tuple<int, int,float>* WeightedEdgeGraph::makeEdgesArray(){
    edgesArray = new std::tuple<int, int, float>[numberOfEdges];
    for (int i =0 ; i<numberOfEdges; i++) {
        edgesArray[i] = edgesVector.at(i);
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

std::vector<std::tuple<int, int,float>> WeightedEdgeGraph::getEdgesVector()const{
    return edgesVector;
}

std::tuple<int, int,float>* WeightedEdgeGraph::getEdgesArray()const{
    return edgesArray;
}

bool WeightedEdgeGraph::adjNodes(int node1, int node2){
    return ( (adjList[node1].find(node2) != adjList[node1].end()) || (adjList[node2].find(node1) != adjList[node2].end())) ;
}