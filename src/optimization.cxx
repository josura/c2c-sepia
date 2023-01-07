#include "Matrix.h"
#include <chrono>
#include <functional>
#include <stdexcept>

template<typename F, typename M>
double run(F f, M const& a, M const& b, M& c)
{
    double time = 0;
    for (int i = 0; i < 10; i++) {
            auto start = std::chrono::high_resolution_clock::now();
            f(a,b,c);
            auto stop = std::chrono::high_resolution_clock::now();

            auto duration = (stop - start);

            time += duration.count();
    }
    return time / 10;
}

template<typename T>
Matrix<T> MatMul(const Matrix<T>& lhs,const Matrix<T>& rhs,Matrix<T>& result){
    //TODO controls over the feasibility of the multiplication
    if(lhs.getCols() == rhs.getRows()){
        for (int r = 0; r < lhs.getRows(); r++) {
                  for (int c = 0; c < rhs.getCols(); c++) {
                          //result(r,c) = 0;   //already initialized at 0
                          for (int i = 0; i < rhs.getRows(); i++)
                                  result(r,c) += lhs.getValue(r,i) * rhs.getValue(i,c);
                  }
          }
    return result;
    } else {
        throw std::invalid_argument("column dimension of lhs is not equal to row dimension of lhs\n");
    }
    
}