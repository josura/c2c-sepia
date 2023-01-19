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
        double* nodeValues;  //arrays of nodeValues
        std::vector<std::unordered_set<int>> adjList;  //adjList as vector of unordered sets
        std::vector<std::string> nameVector;
        std::map<std::string, int> nodeToIndex;

    public:

        //public since they will be accessed a lot
        Matrix<double> adjMatrix;
        std::vector<std::tuple<int, int, double> > edgesVector;

        WeightedEdgeGraph();

        WeightedEdgeGraph(int numNodes);

        WeightedEdgeGraph(std::vector<std::string>& nodeNames);
        WeightedEdgeGraph(std::vector<std::string>& nodeNames,std::vector<double>& nodeVal);

        WeightedEdgeGraph(const Matrix<double>& _adjMatrix);

        ~WeightedEdgeGraph();

        WeightedEdgeGraph* addEdge(int node1, int node2, double weight);
        WeightedEdgeGraph* addEdge(std::string node1name, std::string node2name, double weight);

        double getEdgeWeight(int node1, int node2)const{return adjMatrix.getValue(node1,node2);}
        double getEdgeWeight(std::string node1, std::string node2)const{return adjMatrix.getValue(nodeToIndex.at(node1),nodeToIndex.at(node2));}
        int getIndexFromName(std::string name)const {return nodeToIndex.at(name);}

        WeightedEdgeGraph* addNode(double value);
        WeightedEdgeGraph* addNode(std::string name, double value);


        WeightedEdgeGraph* addNodes(std::vector<double>& values);
        WeightedEdgeGraph* addNodes(std::vector<std::string>& names, std::vector<double>& values);

        WeightedEdgeGraph* setNodeValue(int node, double value);
        WeightedEdgeGraph* setNodeValue(std::string node, double value);

        double getNodeValue(int node)const;
        double getNodeValue(std::string node)const;
        std::vector<double> getNodeValues(std::vector<int> node)const;
        std::vector<double> getNodeValues(std::vector<std::string> node)const;
//TODO controls over nodes and other things

        //optimization functions to make new Matrix, SUGGESTED not using these functions
        Matrix<double> makeMatrix();

        // accessory functions

        int getNumNodes()const ;
        int getNumEdges()const ;
        int degreeOfNode(int node)const;


        std::string getnodeValuesStr()const;

        std::unordered_set<int> getAdjList(int node)const;

        std::string getAdjListStr(int node)const;

        bool adjNodes(int node1, int node2);
        bool adjNodes(std::string node1, std::string node2);

        std::vector<std::tuple<int, int, double>> getEdgesVector()const;



        WeightedEdgeGraph& operator=(const WeightedEdgeGraph& g2);
        //WeightedEdgeGraph& setAdjMatrix(Matrix<double>& mat);

        // optimization methods

        double costFunction(bool* NodeSubset);

        std::vector<int> getSharedAdjacentNodes(std::vector<int>& nodes);

        int getMaxDegree()const;
        double getAverageDegree()const;

        std::map<std::string, int> getNodeToIndexMap()const {return nodeToIndex;}

};


std::ostream& operator<< (std::ostream &out, const WeightedEdgeGraph& data); 