#include <iostream>
#include <vector>
#include "WeightedEdgeGraph.h"
#include "utilities.h"
#include "optimization.cxx"
#include "Matrix.h"

using namespace std;

int main() {

    auto matrix1 = Matrix<double>::createRandom(4,8);
    auto matrix2 = Matrix<double>::createRandom(8,5);
    auto matrixres=Matrix<double>(4,5);
    std::cout<< "runtime for MatMul: " << run(MatMul<double>, matrix1, matrix2, matrixres)<<std::endl;
    // for (int i=0;i<matrixres.getRows() ;i++) {
    //     for (int j=0;j<matrixres.getCols() ;j++) {
    //         cout << matrixres.getValue(i,j) << " ";
    //     }
    //     cout << endl;
    // }
    
    //vector multiplication
    vector<vector<int>> A = {{1, 2}, {3, 4}};
    vector<vector<int>> B = {{5, 6}, {7, 8}};
    vector<vector<int>> result = {{0, 0}, {0, 0}};
    std::cout<< "runtime for MatMulVector(2 x 2 x 2): " << run(matrixMultiplyVector<int>, A, B, result)<<std::endl;
    vector<vector<double>> C = createRandomMatrixVector<double>(4, 8);
    vector<vector<double>> D = createRandomMatrixVector<double>(8, 5);
    vector<vector<double>> result2 = vector<vector<double>>(4,vector<double>(5,0.0));
    std::cout<< "runtime for MatMulVector(4 x 8 x 5): " << run(matrixMultiplyVector<double>, C, D, result2)<<std::endl;
    
    return 0;
}