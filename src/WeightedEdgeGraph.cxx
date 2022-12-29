#pragma once

#include "WeightedEdgeGraph.h"
#include <ostream>
#include <string>
#include <sys/types.h>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <utility>


using uint = unsigned int;

class WeightedEdgeGraph{
    private:
        uint numberOfNodes;
        uint numberOfEdges=0;
        double* nodeWeights;
        std::unordered_set<uint>* adjList;
        std::vector<uint>* adjVector;
        std::vector<std::tuple<uint, uint, float> > edgesVector;
        std::tuple<uint, uint, float>* edgesArray;

    public:
        WeightedEdgeGraph();

        WeightedEdgeGraph(uint numNodes);

        ~WeightedEdgeGraph();

        WeightedEdgeGraph* addEdge(uint node1, uint node2, float weight);

        std::pair<uint, uint>* makeEdgesArray();

        // accessory functions

        uint getNumNodes()const ;
        uint getNumEdges()const ;
        uint degreeOfNode(uint node)const;

        double* getNodeWeights()const;

        std::string getNodeWeightsStr()const;

        std::unordered_set<uint>* getAdjList(uint node)const;

        std::string getAdjListStr(uint node)const;

        bool adjNodes(uint node1, uint node2);

        std::vector<std::pair<uint, uint>> getEdgesVector()const;

        std::pair<uint, uint>* getEdgesArray()const;

        double getNodeWeight(uint node)const;

        // optimization methods

        double costFunction(bool* NodeSubset);

        std::vector<uint> getSharedAdjacentNodes(std::vector<uint>& nodes);

        uint getMaxDegree()const;
        double getAverageDegree()const;

};