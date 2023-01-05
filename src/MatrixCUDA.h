#pragma once

#include "Matrix.h"
#include <iostream>
#include <vector>

template <typename T>
class MatrixCUDA : public Matrix<T> {
    public:
        MatrixCUDA(int, int);
        MatrixCUDA(T**, int, int);
        MatrixCUDA();
        MatrixCUDA(const MatrixCUDA&);
        ~MatrixCUDA();
        MatrixCUDA& operator=(const MatrixCUDA&);

        inline T& operator()(int x, int y) { return _matrix[x][y]; }

        MatrixCUDA& operator+=(const MatrixCUDA&);
        MatrixCUDA& operator-=(const MatrixCUDA&);
        MatrixCUDA& operator*=(const MatrixCUDA&);
        MatrixCUDA& operator*=(T);
        MatrixCUDA& operator/=(T);
        MatrixCUDA operator^(int);  //integer power
        MatrixCUDA operator+(const Matrix<T>&);
        MatrixCUDA operator-(const Matrix<T>&);
        MatrixCUDA operator*(const Matrix<T>&);
        MatrixCUDA operator*(T);
        template<typename U>
        friend Matrix<U> operator*(double, const Matrix<U>&);
        MatrixCUDA operator/(double);
        
        template<typename U>
        friend std::ostream& operator<<(std::ostream&, const Matrix<U>&);
        template<typename U>
        friend std::istream& operator>>(std::istream&, Matrix<U>&);

        void swapRows(int, int);
        MatrixCUDA transpose();

        static MatrixCUDA createIdentity(int);
        static MatrixCUDA solve(MatrixCUDA, MatrixCUDA);
        static MatrixCUDA bandSolve(MatrixCUDA, MatrixCUDA, int);

        // functions on vectors
        static double dotProduct(MatrixCUDA, MatrixCUDA);
        MatrixCUDA& operator*=(const std::vector<T>&);

        // functions for reduction and inverse
        static MatrixCUDA augment(MatrixCUDA, MatrixCUDA);
        MatrixCUDA gaussianElimination();
        MatrixCUDA rowReduceGaussian();
        void readSolutionsFromRREF(std::ostream& os);
        MatrixCUDA inverse();  // to implement

    private:
        int rows_, cols_;
        T **_matrix;

        void allocateMatrixSpace();
        MatrixCUDA expHelper(const MatrixCUDA&, int);
};