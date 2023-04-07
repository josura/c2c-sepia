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