#include "Matrix.h"

template<typename F, typename M>
float run(F f, M const& a, M const& b, M& c);

template<typename T>
Matrix<T> MatMul(const Matrix<T>& lhs,const Matrix<T>& rhs);