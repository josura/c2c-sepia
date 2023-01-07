#include <iostream>
#include "WeightedEdgeGraph.h"
#include "utilities.h"
#include "optimization.cxx"
#include "Matrix.h"

int main() {

    auto matrix1 = Matrix<double>::createRandom(4,8);
    std::cout<<"testing1\n";
    auto matrix2 = Matrix<double>::createRandom(8,5);
    std::cout<<"testing2\n";
    auto matrixres=Matrix<double>(4,5);
    std::cout<<"testing3\n";
    std::cout<< "runtime for MatMul: " << run(MatMul<double>, matrix1, matrix2, matrixres)<<std::endl;
    return 0;
}