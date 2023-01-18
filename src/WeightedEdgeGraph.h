#pragma once

#include "Matrix.h"  //circular dependency :'(
#include <map>
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
        std::vector<std::string> nameVector;
        std::tuple<int, int, double>* edgesArray;
        std::map<std::string, int> nodeToIndex;

    public:

        //public since they will be accessed a lot
        Matrix<double> adjMatrix;
        std::vector<std::tuple<int, int, double> > edgesVector;

        WeightedEdgeGraph();

        WeightedEdgeGraph(int numNodes);

        WeightedEdgeGraph(const Matrix<double>& _adjMatrix);

        ~WeightedEdgeGraph();

        WeightedEdgeGraph* addEdge(int node1, int node2, double weight);
        WeightedEdgeGraph* addEdge(std::string node1name, std::string node2name, double weight);

        WeightedEdgeGraph* addNode(double value);
        WeightedEdgeGraph* addNode(std::string name, double value);


        WeightedEdgeGraph* addNodes(std::vector<double> values);
        WeightedEdgeGraph* addNodes(std::vector<std::string> names, std::vector<double> values);

        WeightedEdgeGraph* setNodeValue(int node, double value);
        WeightedEdgeGraph* setNodeValue(std::string node, double value);
//TODO controls over nodes and other things

        //optimization functions to make EdgesArray and new Matrix, SUGGESTED not using these functions
        std::tuple<int, int, double>* makeEdgesArray(); //edge =(node, node, weight)
        Matrix<double> makeMatrix();

        // accessory functions

        int getNumNodes()const ;
        int getNumEdges()const ;
        int degreeOfNode(int node)const;

        double* getNodeWeights()const;

        std::string getNodeWeightsStr()const;

        std::unordered_set<int>* getAdjList(int node)const;

        std::string getAdjListStr(int node)const;

        bool adjNodes(int node1, int node2);
        bool adjNodes(std::string node1, std::string node2);

        std::vector<std::tuple<int, int, double>> getEdgesVector()const;

        std::tuple<int, int, double>* getEdgesArray()const;

        double getNodeWeight(int node)const;


        WeightedEdgeGraph& operator=(const WeightedEdgeGraph& g2);
        //WeightedEdgeGraph& setAdjMatrix(Matrix<double>& mat);

        // optimization methods

        double costFunction(bool* NodeSubset);

        std::vector<int> getSharedAdjacentNodes(std::vector<int>& nodes);

        int getMaxDegree()const;
        double getAverageDegree()const;

};


std::ostream& operator<< (std::ostream &out, const WeightedEdgeGraph& data); 