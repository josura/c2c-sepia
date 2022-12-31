#include "utilities.h"
#include "WeightedEdgeGraph.h"
#include <cstddef>
#include <iostream>
#include <random>
#include <tuple>
#include <unordered_set>
#include <vector>

using namespace std;

std::ostream& operator<< (std::ostream &out, WeightedEdgeGraph const& data) {
            out << data.getNumNodes() << " " << data.getNumEdges() <<std::endl;
            string nodeweights = data.getNodeWeightsStr();
            out << nodeweights << std::endl;
            out << "Adj Lists" << std::endl;
            for(int i = 0; i<data.getNumNodes(); i++){
                out << "node " << i << " :" << data.getAdjListStr(i) << std::endl;
            }
            out << "Edges vector: {";
            for(auto it = data.getEdgesVector()->cbegin();it!=data.getEdgesVector()->cend();it++){
                out << "(" << get<0>(*it)<< ","<< get<1>(*it)<< ","<< get<2>(*it)<< "," << ")," ;
            }
            out << "}"<< std::endl;
            return out;
        }
