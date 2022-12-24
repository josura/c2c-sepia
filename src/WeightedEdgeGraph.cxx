#pragma once

#include <ostream>
#include <string>
#include <sys/types.h>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <utility>

class WeightedEdgeGraph{
    private:
        uint numberOfNodes;
        uint numberOfEdges=0;
        double* nodeWeights;
        std::unordered_set<uint>* adjList;
        std::vector<uint>* adjVector;
        std::vector<std::pair<uint, uint> > edgesVector;
        std::pair<uint, uint>* edgesArray;

    public:
        WeightedEdgeGraph();

        WeightedEdgeGraph(uint numNodes, double* nodeWeights);

        ~WeightedEdgeGraph();

        WeightedEdgeGraph* addEdge(uint node1, uint node2);

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