#include <armadillo>
#include <vector>

template<typename T>
std::vector<std::vector<T>> armaMatrixToVector(arma::Mat<T> matr);

template<typename T>
std::vector<T> armaColumnToVector(arma::Col<T> matr);

template<typename T>
std::vector<T> armaRowToVector(arma::Row<T> matr);