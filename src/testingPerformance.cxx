#include <iostream>
#include "WeightedEdgeGraph.h"
#include "utilities.h"
#include "optimization.cxx"
#include "Matrix.h"

int main() {

    auto matrix1 = Matrix<double>::createRandom(4,8);
    
    auto matrix2 = Matrix<double>::createRandom(8,5);
    auto matrixres=Matrix<double>(4,5);
    run(MatMul<double>, matrix1, matrix2, matrixres);
    return 0;
}