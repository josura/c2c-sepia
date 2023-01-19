#include <iostream>
#include "WeightedEdgeGraph.h"
#include "utilities.h"

int main() {
    WeightedEdgeGraph* g1 = new WeightedEdgeGraph(4);
    WeightedEdgeGraph* g2 = new WeightedEdgeGraph();
    WeightedEdgeGraph gtest = *g1;
    std::cout << " testing "<< g2->getNumEdges()<<std::endl;
    std::cout << *g2 ;
    g1->addEdge(1,2,0.1)->addEdge(1,3,0.4);
    *g2 = *g1;
    WeightedEdgeGraph graphtest = *g1;
    std::cout << " testing "<< g1->getNumEdges()<<std::endl;
    std::cout << *g1 ;
    std::cout << " testing "<< graphtest.getNumEdges()<<std::endl;
    std::cout << graphtest ;
    std::cout << " testing "<< g2->getNumEdges()<<std::endl;
    std::cout << *g2 ;
    return 0;
}