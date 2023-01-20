#pragma once

#include "Matrix.h"
#include "WeightedEdgeGraph.h"
#include <map>
#include <string>
#include <tuple>
#include <vector>

class Computation{
    private:
        std::vector<double> input,output;
        WeightedEdgeGraph* metapathway;
        WeightedEdgeGraph* augmentedMetapathway;
        std::map<std::string, int> geneMapToNode;
        std::vector<std::string> cellTypes;
        std::string localCellType;
    public:
        Computation();
        Computation(std::string _thisCellType, const std::vector<double>& _input);   // useless???

        /*
        Computation without knowledge of the other cell types, this part can be seen as the classical algorithm without additional computation for message passing
        between cells, only intra-cell propagation
        @param std::string _thisCellType: the celltype of this computation, this information will be used as the unique name for the Agent
        @param const std::vector<double>& _input: input vector of the nodes values, initially the one passed in the input
        @param const Matrix<double>& _W: the adjacency matrix along the values of every edge in the graph that it represents
        @param const std::vector<std::string>& metapathwayNames: the graph nodes names, in order defined by the adjacency matrix
        */
        Computation(std::string _thisCellType,const std::vector<double>& _input, const Matrix<double>& _W, const std::vector<std::string>& metapathwayNames);
        
        /*
        Computation with the knowledge of the other cell types, intra-cell and inter-cell propagation
        @param std::string _thisCellType: the celltype of this computation, this information will be used as the unique name for the Agent and to filter the other celltypes
        @param const std::vector<double>& _input: input vector of the nodes values, initially the one passed in the input
        @param const Matrix<double>& _W: the adjacency matrix along the values of every edge in the graph that it represents
        @param const std::vector<std::string>& metapathwayNames: the graph nodes names, in order defined by the adjacency matrix
        @param const std::vector<std::string>& _cellTypes: the celltypes other than this celltype, the other agents in the network
        */
        Computation(std::string _thisCellType,const std::vector<double>& _input, const Matrix<double>& _W, const std::vector<std::string>& metapathwayNames, const std::vector<std::string>& _cellTypes);
        void augmentMetapathway(std::vector<std::string>&,std::vector<std::tuple<std::string,std::string,double>>&);
        std::vector<double> computePerturbation();
        std::vector<double> computeAugmentedPerturbation(); //taking into account virtual nodes in the augmented metapathway

        
};