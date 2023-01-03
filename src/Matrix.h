#pragma once

#include <iostream>
#include <vector>

template <typename T>
class Matrix {
    public:
        Matrix(int, int);
        Matrix(T**, int, int);
        Matrix();
        Matrix(const Matrix&);
        ~Matrix();
        Matrix& operator=(const Matrix&);

        inline T& operator()(int x, int y) { return _matrix[x][y]; }

        Matrix& operator+=(const Matrix&);
        Matrix& operator-=(const Matrix&);
        Matrix& operator*=(const Matrix&);
        Matrix& operator*=(T);
        Matrix& operator/=(T);
        Matrix operator^(int);  //integer power
        Matrix operator+(const Matrix<T>&);
        Matrix operator-(const Matrix<T>&);
        Matrix operator*(const Matrix<T>&);
        Matrix operator*(T);
        template<typename U>
        friend Matrix<U> operator*(double, const Matrix<U>&);
        Matrix operator/(double);
        
        template<typename U>
        friend std::ostream& operator<<(std::ostream&, const Matrix<U>&);
        template<typename U>
        friend std::istream& operator>>(std::istream&, Matrix<U>&);

        void swapRows(int, int);
        Matrix transpose();

        static Matrix createIdentity(int);
        static Matrix solve(Matrix, Matrix);
        static Matrix bandSolve(Matrix, Matrix, int);

        // functions on vectors
        static double dotProduct(Matrix, Matrix);
        Matrix& operator*=(const std::vector<T>&);

        // functions for reduction and inverse
        static Matrix augment(Matrix, Matrix);
        Matrix gaussianElimination();
        Matrix rowReduceGaussian();
        void readSolutionsFromRREF(std::ostream& os);
        Matrix inverse();  // to implement

    private:
        int rows_, cols_;
        T **_matrix;

        void allocateMatrixSpace();
        Matrix expHelper(const Matrix&, int);
};