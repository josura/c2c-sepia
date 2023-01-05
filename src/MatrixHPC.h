#pragma once

#include "Matrix.h"
#include <iostream>
#include <vector>

template <typename T>
class MatrixHPC : public Matrix<T> {
    public:
        MatrixHPC(int, int);
        MatrixHPC(T**, int, int);
        MatrixHPC();
        MatrixHPC(const  MatrixHPC&);
        ~ MatrixHPC();
        MatrixHPC& operator=(const  MatrixHPC&);

        inline T& operator()(int x, int y) { return _matrix[x][y]; }

        MatrixHPC& operator+=(const  MatrixHPC&);
        MatrixHPC& operator-=(const  MatrixHPC&);
        MatrixHPC& operator*=(const  MatrixHPC&);
        MatrixHPC& operator*=(T);
        MatrixHPC& operator/=(T);
        MatrixHPC operator^(int);  //integer power
        MatrixHPC operator+(const Matrix<T>&);
        MatrixHPC operator-(const Matrix<T>&);
        MatrixHPC operator*(const Matrix<T>&);
        MatrixHPC operator*(T);
        template<typename U>
        friend Matrix<U> operator*(double, const Matrix<U>&);
        MatrixHPC operator/(double);
        
        template<typename U>
        friend std::ostream& operator<<(std::ostream&, const Matrix<U>&);
        template<typename U>
        friend std::istream& operator>>(std::istream&, Matrix<U>&);

        void swapRows(int, int);
        MatrixHPC transpose();

        static MatrixHPC createIdentity(int);
        static MatrixHPC solve(MatrixHPC, MatrixHPC);
        static MatrixHPC bandSolve(MatrixHPC, MatrixHPC, int);

        // functions on vectors
        static double dotProduct(MatrixHPC, MatrixHPC);
        MatrixHPC& operator*=(const std::vector<T>&);

        // functions for reduction and inverse
        static MatrixHPC augment(MatrixHPC, MatrixHPC);
        MatrixHPC gaussianElimination();
        MatrixHPC rowReduceGaussian();
        void readSolutionsFromRREF(std::ostream& os);
        MatrixHPC inverse();  // to implement

    private:
        int rows_, cols_;
        T **_matrix;

        void allocateMatrixSpace();
        MatrixHPC expHelper(const  MatrixHPC&, int);
};