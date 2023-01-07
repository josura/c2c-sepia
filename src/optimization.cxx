#include "Matrix.h"
#include <chrono>
#include <functional>

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
    for (int r = 0; r < lhs.dim(); r++) {
                  for (int c = 0; c < lhs.dim(); c++) {
                          result(r,c) = 0;
                          for (int i = 0; i < lhs.dim(); i++)
                                  result(r,c) += lhs(r,i) * rhs(i,c);
                  }
          }
    return result;
}