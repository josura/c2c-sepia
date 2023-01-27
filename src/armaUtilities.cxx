#include "armaUtilities.h"
#include <vector>

template<typename T>
std::vector<std::vector<T>> armaMatrixToVector(arma::Mat<T> matr){
    std::vector<std::vector<T>> ret(matr.n_rows);
    for (size_t i = 0; i < matr.n_rows; ++i) {
        ret[i] = arma::conv_to< std::vector<T> >::from(matr.row(i));
    };
    return ret;
}

template<typename T>
std::vector<T> armaColumnToVector(arma::Col<T> matr){
    return arma::conv_to< std::vector<T> >::from(matr);
}

template<typename T>
std::vector<T> armaRowToVector(arma::Row<T> matr){
    return arma::conv_to< std::vector<T> >::from(matr);
}