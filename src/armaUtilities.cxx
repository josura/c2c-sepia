#include<iostream>
#include "armaUtilities.h"
#include <vector>

void print_mat(arma::mat my_matrix){
    
    uint cols = my_matrix.n_cols;
    uint rows = my_matrix.n_rows;
    
    std::cout << "--------\n";
    for(uint rX = 0; rX < rows; rX++) {
        for(uint cX = 0; cX < cols; cX++) {
            std::cout << my_matrix(rX, cX) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "--------\n";
}

template<typename T>
arma::Mat<T> normalizeColumns(arma::Mat<T> matr){
    arma::Mat<T> normalizedMatr(matr.n_rows, matr.n_cols);
    for (uint cX = 0; cX < matr.n_cols; cX++) {
        arma::Col<T> col = matr.col(cX);
        normalizedMatr.col(cX) = col / (arma::norm(col) + 1e-10);
    }
    return normalizedMatr;
}

template arma::Mat<double> normalizeColumns(arma::Mat<double> matr);


template<typename T>
arma::Mat<T> normalizeRows(arma::Mat<T> matr){
    arma::Mat<T> normalizedMatr(matr.n_rows, matr.n_cols);
    for (uint rX = 0; rX < matr.n_rows; rX++) {
        arma::Row<T> row = matr.row(rX);
        normalizedMatr.row(rX) = row / (arma::norm(row) + 1e-10);
    }
    return normalizedMatr;
}

template arma::Mat<double> normalizeRows(arma::Mat<double> matr);


template<typename T>
arma::Mat<T> normalize1Columns(arma::Mat<T> matr){
    arma::Mat<T> normalizedMatr(matr.n_rows, matr.n_cols);
    for (uint cX = 0; cX < matr.n_cols; cX++) {
        arma::Col<T> col = matr.col(cX);
        normalizedMatr.col(cX) = col / (arma::sum(arma::abs(col)) + 1e-10);
    }
    return normalizedMatr;
}

template arma::Mat<double> normalize1Columns(arma::Mat<double> matr);


template<typename T>
arma::Mat<T> normalize1Rows(arma::Mat<T> matr){
    arma::Mat<T> normalizedMatr(matr.n_rows, matr.n_cols);
    for (uint rX = 0; rX < matr.n_rows; rX++) {
        arma::Row<T> row = matr.row(rX);
        normalizedMatr.row(rX) = row / (arma::sum(arma::abs(row)) + 1e-10);
    }
    return normalizedMatr;
}

template arma::Mat<double> normalize1Rows(arma::Mat<double> matr);