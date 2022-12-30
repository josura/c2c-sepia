#include <iostream>
#include "WeightedEdgeGraph.h"

int main() {
    std::cout << "c2c-sepia init";
    WeightedEdgeGraph* g1 = new WeightedEdgeGraph(4);
    g1->addEdge(1,2,0.1)->addEdge(1,3,0.4);
    std::cout << "testing \n"<< g1->getNumEdges();
    return 0;
}