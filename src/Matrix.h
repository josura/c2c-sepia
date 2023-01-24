#pragma once

#include <iostream>
#include <iterator>
#include <vector>
#include "utilities.h"

template <typename T>
class Matrix {
    public:
        Matrix(int, int);
        Matrix(T**, int, int);
        Matrix();
        Matrix(const Matrix&);
        Matrix(const T*);
        Matrix(const std::vector<T>&);
        ~Matrix();
        Matrix& operator=(const Matrix&);

        inline T& operator()(int x, int y) { return _matrix[x][y]; }
        inline T& getValue(int x, int y)const{ return _matrix[x][y]; }

        Matrix& operator+=(const Matrix&);
        Matrix& operator-=(const Matrix&);
        Matrix& operator*=(const Matrix&);
        Matrix& operator*=(T);
        Matrix& operator/=(T);
        Matrix operator^(int);  //integer power
        Matrix operator+(const Matrix<T>&)const;
        Matrix operator-(const Matrix<T>&)const;
        Matrix operator*(const Matrix<T>&)const;
        Matrix operator*(const std::vector<T>&)const;
        Matrix operator*(T)const;
        template<typename U>
        friend Matrix<U> operator*(double, const Matrix<U>&);
        template<typename U>
        friend Matrix<U> operator*(U*,const Matrix<U>&);  //vector multiplication leftwise
        template<typename U>
        friend Matrix<U> operator*(std::vector<U>&,const Matrix<U>&);
        Matrix operator/(double);
        
        template<typename U>
        friend std::ostream& operator<<(std::ostream&, const Matrix<U>&);
        template<typename U>
        friend std::istream& operator>>(std::istream&, Matrix<U>&);

        void swapRows(int, int);
        Matrix transpose()const;

        static Matrix createIdentity(int);
        static Matrix createRandom(int,int);
        static Matrix solve(Matrix, Matrix);
        static Matrix bandSolve(Matrix, Matrix, int);

        // functions on vectors
        static double dotProduct(Matrix, Matrix);
        static Matrix getMinor(const Matrix<T>&,int, int,int);
        static T determinant(const Matrix<T>& A);
        Matrix& operator*=(const std::vector<T>&);

        // functions for reduction and inverse
        Matrix concatenateRight(const Matrix&)const;
        T determinant()const;
        Matrix gaussianElimination();
        Matrix rowReduceGaussian();
        Matrix inverse();  // to implement

        //get functions
        int getRows()const{return rows_;}
        int getCols()const{return cols_;}
        bool isVector()const{return (rows_ > 1 && cols_ == 1); }
        std::vector<T> asVector()const; 

        //functions to add rows and columns while mantaining the original data in the upperleft corner(these functions are bad, better use a vector when trying to work with dynamically instantiated data)
        // also these functions create a copy and do not work on the original
        Matrix copyAndAddRowsCols(int additionalRows, int additionalCols) const;

    private:
        int rows_, cols_;
        T **_matrix;

        void allocateMatrixSpace();
        Matrix expHelper(const Matrix&, int);
};