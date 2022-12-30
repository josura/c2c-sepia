#pragma once

#include <ostream>
#include <string>
#include <sys/types.h>
#include <tuple>
/*std::get<i>( tpl) to get the value or set it in a tuple*/
#include <unordered_set>
#include <iostream>
#include <vector>
#include <utility>

class WeightedEdgeGraph{
    private:
        int numberOfNodes;
        int numberOfEdges=0;
        double* nodeWeights;
        std::unordered_set<int>* adjList;
        std::vector<int>* adjVector;
        std::vector<std::tuple<int, int, float> > edgesVector;
        std::tuple<int, int, float>* edgesArray;

    public:
        WeightedEdgeGraph();

        WeightedEdgeGraph(int numNodes);

        ~WeightedEdgeGraph();

        WeightedEdgeGraph* addEdge(int node1, int node2, float weight);

        std::tuple<int, int, float>* makeEdgesArray(); //edge =(node, node, weight)

        // accessory functions

        int getNumNodes()const ;
        int getNumEdges()const ;
        int degreeOfNode(int node)const;

        double* getNodeWeights()const;

        std::string getNodeWeightsStr()const;

        std::unordered_set<int>* getAdjList(int node)const;

        std::string getAdjListStr(int node)const;

        bool adjNodes(int node1, int node2);

        std::vector<std::tuple<int, int, float>> getEdgesVector()const;

        std::tuple<int, int, float>* getEdgesArray()const;

        double getNodeWeight(int node)const;

        // optimization methods

        double costFunction(bool* NodeSubset);

        std::vector<int> getSharedAdjacentNodes(std::vector<int>& nodes);

        int getMaxDegree()const;
        double getAverageDegree()const;

};