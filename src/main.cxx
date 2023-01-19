#include <iostream>
#include "WeightedEdgeGraph.h"
#include "utilities.h"

int main() {
    WeightedEdgeGraph* g1 = new WeightedEdgeGraph(4);
    WeightedEdgeGraph* g2 = new WeightedEdgeGraph();
    std::cout << " testing "<< g2->getNumEdges()<<std::endl;
    std::cout << *g2 ;
    g1->addEdge(1,2,0.1)->addEdge(1,3,0.4);
    *g2 = *g1;
    WeightedEdgeGraph graphtest = *g1;
    std::cout << " testing "<< g1->getNumEdges()<<std::endl;
    std::cout << *g1 ;
    std::cout << " testing "<< g2->getNumEdges()<<std::endl;
    std::cout << *g2 ;

    std::vector<std::string> nodeNames{"node1","node2","node3","node4","node5"};
    std::vector<double> nodeValues{0.3,4.1,3.8,8.2,9.5};
    WeightedEdgeGraph g4 = WeightedEdgeGraph(nodeNames);
    WeightedEdgeGraph g5 = WeightedEdgeGraph(nodeNames,nodeValues);
    WeightedEdgeGraph g6;
    g6.assign(g5);

    g5.addEdge("node1","node3",2.3)->addEdge("node3","node2",0.3);

    std::cout << " testing "<< g4.getNumEdges()<<std::endl;
    std::cout << g4 ;
    std::cout << " testing "<< g5.getNumEdges()<<std::endl;
    std::cout << g5 ;
    //std::cout << " testing "<< g6.getNumEdges()<<std::endl;
    //std::cout << g6 ;
    return 0;
}