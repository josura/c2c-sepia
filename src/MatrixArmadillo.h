#pragma once

#include <iostream>
#include <iterator>
#include <vector>
#include "utilities.h"

template <typename T>
class MatrixArmadillo {
    public:
        MatrixArmadillo(int, int);
        MatrixArmadillo(T**, int, int);
        MatrixArmadillo();
        MatrixArmadillo(const MatrixArmadillo&);
        MatrixArmadillo(const T*);
        MatrixArmadillo(const std::vector<T>&);
        ~MatrixArmadillo();
        MatrixArmadillo& operator=(const MatrixArmadillo&);

        inline T& operator()(int x, int y) { return _matrix[x][y]; }
        inline T& getValue(int x, int y)const{ return _matrix[x][y]; }

        MatrixArmadillo& operator+=(const MatrixArmadillo&);
        MatrixArmadillo& operator-=(const MatrixArmadillo&);
        MatrixArmadillo& operator*=(const MatrixArmadillo&);
        MatrixArmadillo& operator*=(T);
        MatrixArmadillo& operator/=(T);
        MatrixArmadillo operator^(int);  //integer power
        MatrixArmadillo operator+(const MatrixArmadillo<T>&)const;
        MatrixArmadillo operator-(const MatrixArmadillo<T>&)const;
        MatrixArmadillo operator*(const MatrixArmadillo<T>&)const;
        MatrixArmadillo operator*(const std::vector<T>&)const;
        MatrixArmadillo operator*(T)const;
        template<typename U>
        friend MatrixArmadillo<U> operator*(double, const MatrixArmadillo<U>&);
        template<typename U>
        friend MatrixArmadillo<U> operator*(U*,const MatrixArmadillo<U>&);  //vector multiplication leftwise
        template<typename U>
        friend MatrixArmadillo<U> operator*(std::vector<U>&,const MatrixArmadillo<U>&);
        MatrixArmadillo operator/(double);
        
        template<typename U>
        friend std::ostream& operator<<(std::ostream&, const MatrixArmadillo<U>&);
        template<typename U>
        friend std::istream& operator>>(std::istream&, MatrixArmadillo<U>&);

        void swapRows(int, int);
        MatrixArmadillo transpose()const;

        static MatrixArmadillo createIdentity(int);
        static MatrixArmadillo createRandom(int,int);
        static MatrixArmadillo solve(MatrixArmadillo, MatrixArmadillo);
        static MatrixArmadillo bandSolve(MatrixArmadillo, MatrixArmadillo, int);

        // functions on vectors
        static double dotProduct(MatrixArmadillo, MatrixArmadillo);
        static MatrixArmadillo getMinor(const MatrixArmadillo<T>&,int, int,int);
        static T determinant(const MatrixArmadillo<T>& A);
        MatrixArmadillo& operator*=(const std::vector<T>&);

        // functions for reduction and inverse
        MatrixArmadillo concatenateRight(const MatrixArmadillo&)const;
        T determinant()const;
        MatrixArmadillo gaussianElimination();
        MatrixArmadillo rowReduceGaussian();
        MatrixArmadillo inverse();  // to implement

        //get functions
        int getRows()const{return rows_;}
        int getCols()const{return cols_;}
        bool isVector()const{return (rows_ > 1 && cols_ == 1); }
        std::vector<T> asVector()const; 

        //functions to add rows and columns while mantaining the original data in the upperleft corner(these functions are bad, better use a vector when trying to work with dynamically instantiated data)
        // also these functions create a copy and do not work on the original
        MatrixArmadillo copyAndAddRowsCols(int additionalRows, int additionalCols) const;

    private:
        int rows_, cols_;
        T **_matrix;

        void allocateMatrixSpace();
        MatrixArmadillo expHelper(const MatrixArmadillo&, int);
};