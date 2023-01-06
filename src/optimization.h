#include "Matrix.h"

template<typename T, typename M>
float run(T f, M const& a, M const& b, M& c);

template<typename T>
Matrix<T> MatMul(const Matrix<T>& lhs,const Matrix<T>& rhs);