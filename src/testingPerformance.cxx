#include <iostream>
#include "WeightedEdgeGraph.h"
#include "utilities.h"
#include "optimization.h"
#include "Matrix.h"
int main() {
    WeightedEdgeGraph* g1 = new WeightedEdgeGraph(4);
    WeightedEdgeGraph* g2 = new WeightedEdgeGraph();
    std::cout << " testing "<< g2->getNumEdges()<<std::endl;
    std::cout << *g2 ;
    g1->addEdge(1,2,0.1)->addEdge(1,3,0.4);
    *g2 = *g1;
    std::cout << " testing "<< g1->getNumEdges()<<std::endl;
    std::cout << *g1 ;
    std::cout << " testing "<< g2->getNumEdges()<<std::endl;
    std::cout << *g2 ;

    auto matrix1 = Matrix<double>::createRandom(4,8);
    
    auto matrix2 = Matrix<double>::createRandom(8,5);
    auto matrixres=Matrix<double>(4,5);
    run(MatMul<double>, matrix1, matrix2, matrixres);
    return 0;
}