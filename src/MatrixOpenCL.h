#pragma once

#include "Matrix.h"
#include <iostream>
#include <vector>

template <typename T>
class MatrixOpenCL : public Matrix<T> {
    public:
        MatrixOpenCL(int, int);
        MatrixOpenCL(T**, int, int);
         MatrixOpenCL();
         MatrixOpenCL(const  MatrixOpenCL&);
        ~ MatrixOpenCL();
        MatrixOpenCL& operator=(const  MatrixOpenCL&);

        inline T& operator()(int x, int y) { return _matrix[x][y]; }

        MatrixOpenCL& operator+=(const  MatrixOpenCL&);
        MatrixOpenCL& operator-=(const  MatrixOpenCL&);
        MatrixOpenCL& operator*=(const  MatrixOpenCL&);
        MatrixOpenCL& operator*=(T);
        MatrixOpenCL& operator/=(T);
        MatrixOpenCL operator^(int);  //integer power
        MatrixOpenCL operator+(const Matrix<T>&);
        MatrixOpenCL operator-(const Matrix<T>&);
        MatrixOpenCL operator*(const Matrix<T>&);
        MatrixOpenCL operator*(T);
        template<typename U>
        friend Matrix<U> operator*(double, const Matrix<U>&);
        MatrixOpenCL operator/(double);
        
        template<typename U>
        friend std::ostream& operator<<(std::ostream&, const Matrix<U>&);
        template<typename U>
        friend std::istream& operator>>(std::istream&, Matrix<U>&);

        void swapRows(int, int);
        MatrixOpenCL transpose();

        static MatrixOpenCL createIdentity(int);
        static MatrixOpenCL solve(MatrixOpenCL, MatrixOpenCL);
        static MatrixOpenCL bandSolve(MatrixOpenCL, MatrixOpenCL, int);

        // functions on vectors
        static double dotProduct(MatrixOpenCL, MatrixOpenCL);
         MatrixOpenCL& operator*=(const std::vector<T>&);

        // functions for reduction and inverse
        static MatrixOpenCL augment(MatrixOpenCL, MatrixOpenCL);
        MatrixOpenCL gaussianElimination();
        MatrixOpenCL rowReduceGaussian();
        void readSolutionsFromRREF(std::ostream& os);
        MatrixOpenCL inverse();  // to implement

    private:
        int rows_, cols_;
        T **_matrix;

        void allocateMatrixSpace();
        MatrixOpenCL expHelper(const  MatrixOpenCL&, int);
};