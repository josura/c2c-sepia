#include <armadillo>
#include <vector>

template<typename T>
std::vector<std::vector<T>> armaMatrixToVector(arma::Mat<T> matr){return arma::conv_to< std::vector<T> >::from(arma::vectorise( matr));}

template<typename T>
std::vector<T> armaColumnToVector(arma::Col<T> matr){return arma::conv_to< std::vector<T> >::from(matr);}

template<typename T>
std::vector<T> armaRowToVector(arma::Row<T> matr){return arma::conv_to< std::vector<T> >::from(matr);}

template<typename T>
arma::Col<T> vectorToArmaColumn(std::vector<T> vec){return arma::Col<T>(vec);}

void print_mat(arma::mat my_matrix);