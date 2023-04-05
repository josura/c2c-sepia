#pragma once

#include "Matrix.h"  //circular dependency :'(
#include <map>
#include <ostream>
#include <stdexcept>
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
        int numberOfNodes=0;
        int numberOfEdges=0;
        double* nodeValues=nullptr;  //arrays of nodeValues
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

        double getEdgeWeight(int node1, int node2)const{
            if(node1 >= 0 && node2 >= 0)
                return adjMatrix.getValue(node1,node2);
            else throw std::out_of_range("WeightedEdgeGraph::getEdgeWeight: one of the nodes is out of range(index)");
        }
        double getEdgeWeight(std::string node1, std::string node2)const{
            if(getIndexFromName(node1) >= 0 && getIndexFromName(node2) >= 0)
                return adjMatrix.getValue(nodeToIndex.at(node1),nodeToIndex.at(node2));
            else throw std::out_of_range("WeightedEdgeGraph::getEdgeWeight: one of the nodes is out of range(string)");
        }

        
        //return -1 if the name is not found, so be careful
        int getIndexFromName(std::string name)const {
            try {
                return nodeToIndex.at(name);
            }
            catch (const std::out_of_range&) {
                return -1;
                // TODO: Deal with the missing element.
            }
        }

        WeightedEdgeGraph* addNode(double value=0);
        WeightedEdgeGraph* addNode(std::string name, double value=0);

        /*
        add n nodes, where n is the size of the vector passed as input, and assign to each node its corrispondent value
        */
        WeightedEdgeGraph* addNodes(const std::vector<double>& values);
        /*
        add n nodes with names, where n is the size of the vector passed as input, and assign to each node its corrispondent value
        */
        WeightedEdgeGraph* addNodes(const std::vector<std::string>& names, const std::vector<double>& values=std::vector<double>());

        /*
        add n nodes, where n is the size of the vector passed as input, and assign to each node its corrispondent value
        and copy the graph into a new graph in dynamic memory
        */
        WeightedEdgeGraph* addNodesAndCopyNew(const std::vector<double>& values);
        
        /*
        add n nodes with names, where n is the size of the vector passed as input, and assign to each node its corrispondent value
        and copy the graph into a new graph in dynamic memory
        */
        WeightedEdgeGraph* addNodesAndCopyNew(const std::vector<std::string>& names, const std::vector<double>& values=std::vector<double>());


        /*
        set node values by index
        */
        WeightedEdgeGraph* setNodeValue(int node, double value);
        
        /*
        set node values by name
        */
        WeightedEdgeGraph* setNodeValue(std::string node, double value);
        //WeightedEdgeGraph* setNodesValues(std::vector< double> value,const std::vector<std::string> node = std::vector<std::string>());

        /*
        change node name
        */
        WeightedEdgeGraph* setNodeName(std::string nodenameTarget, std::string nodenameSet);
        /*
        If provided with one parameter, controls the vector size and sets the node names if they are of the same size.
        If provided with two parameters, changes the nodes in nodenameTargets with the values in nodenameSets 
        */
        WeightedEdgeGraph* setNodesNames(const std::vector<std::string>& nodenameSets, const std::vector<std::string>& nodenameTargets = std::vector<std::string>());

        double getNodeValue(int node)const;
        double getNodeValue(std::string node)const;
        std::vector<double> getNodeValues(const std::vector<int>& node)const;
        std::vector<double> getNodeValues(const std::vector<std::string>& node=std::vector<std::string>())const;
//TODO controls over nodes and other things

        //optimization functions to make new Matrix, SUGGESTED not using these functions
        Matrix<double> makeMatrix();

        // accessory functions

        int getNumNodes()const ;
        int getNumEdges()const ;
        std::vector<std::string> getNodeNames()const{ return nameVector;}
        int degreeOfNode(int node)const;


        std::string getnodeValuesStr()const;

        std::unordered_set<int> getAdjList(int node)const;
        std::unordered_set<int> getAdjList(std::string node)const;

        std::string getAdjListStr(int node)const;
        std::string getAdjListStr(std::string node)const;


        //return true if the nodes are adjacent (by index)
        bool adjNodes(int node1, int node2);
        //return true if the nodes are adjacent (by name)
        bool adjNodes(std::string node1, std::string node2);
        //return true if the nodes are connected (by index)
        bool connectedNodes(int node1, int node2);
        //return true if the nodes are connected (by name)
        bool connectedNodes(std::string node1, std::string node2);

        std::vector<std::tuple<int, int, double>> getEdgesVector()const;



        WeightedEdgeGraph& operator=(const WeightedEdgeGraph& g2);
        //WeightedEdgeGraph& operator=(const WeightedEdgeGraph g2);
        //WeightedEdgeGraph& setAdjMatrix(Matrix<double>& mat);
        void assign(const WeightedEdgeGraph& g2);
        WeightedEdgeGraph* copyNew()const;

        // optimization methods

        double costFunction(bool* NodeSubset);

        std::vector<int> getSharedAdjacentNodes(std::vector<int>& nodes);

        int getMaxDegree()const;
        double getAverageDegree()const;

        std::map<std::string, int> getNodeToIndexMap()const {return nodeToIndex;}
        void print()const;

};


std::ostream& operator<< (std::ostream &out, const WeightedEdgeGraph& data); 