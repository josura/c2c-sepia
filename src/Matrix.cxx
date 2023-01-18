#include "Matrix.h"
#include "gtest/gtest.h"


// helper functions
template<typename T>
void Matrix<T>::allocateMatrixSpace()
{
    _matrix = new T*[rows_];
    for (int i = 0; i < rows_; i++) {
        _matrix[i] = new T[cols_];
    }
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            _matrix[i][j] = 0;
        }    
    }
}

template void Matrix<double>::allocateMatrixSpace();


//constructors and destructors
template<typename T>
Matrix<T>::Matrix(int rows, int cols) : rows_(rows), cols_(cols)
{
    allocateMatrixSpace();
}

template Matrix<double>::Matrix(int rows, int cols);

template<typename T>
Matrix<T>::Matrix(T** a, int rows, int cols) : rows_(rows), cols_(cols)
{
    allocateMatrixSpace();
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            _matrix[i][j] = a[i][j];
        }
    }
}

template Matrix<double>::Matrix(double** a, int rows, int cols);

template<typename T>
Matrix<T>::Matrix() : rows_(1), cols_(1)
{
    allocateMatrixSpace();
    _matrix[0][0] = 0;
}

template Matrix<double>::Matrix();

template<typename T>
Matrix<T>::Matrix(const std::vector<T>& vec):rows_(vec.size()),cols_(1){
    allocateMatrixSpace();
    for (int i = 0; i < SizeToInt(vec.size()); i++) {
        _matrix[i][0] = vec[i];
    }
}

template Matrix<double>::Matrix(const std::vector<double>& vec);


template<typename T>
Matrix<T>::~Matrix()
{
    for (int i = 0; i < rows_; ++i) {
        delete[] _matrix[i];
    }
    delete[] _matrix;
}

template Matrix<double>::~Matrix();

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& m) : rows_(m.rows_), cols_(m.cols_)
{
    allocateMatrixSpace();
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            _matrix[i][j] = m._matrix[i][j];
        }
    }
}

template Matrix<double>::Matrix(const Matrix<double>& m);

//static methods
template<typename T>
Matrix<T> Matrix<T>::createRandom(int rows,int cols){
    Matrix<T> retMat=Matrix<T>(rows,cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            setRandom(retMat(i,j));
        }
    }
    return retMat;
}

template Matrix<double> Matrix<double>::createRandom(int rows,int cols);

template<typename T>
Matrix<T> Matrix<T>::createIdentity(int size)
{
    Matrix<T> temp(size, size);
    for (int i = 0; i < temp.rows_; ++i) {
        for (int j = 0; j < temp.cols_; ++j) {
            if (i == j) {
                temp(i,j) = 1;
            } else {
                temp(i,j) = 0;
            }
        }
    }
    return temp;
}
template Matrix<double> Matrix<double>::createIdentity(int size);

//operators
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
    if (this == &m) {
        return *this;
    }

    if (rows_ != m.rows_ || cols_ != m.cols_) {
        for (int i = 0; i < rows_; i++) {
            delete[] _matrix[i];
        }
        delete[] _matrix;

        rows_ = m.rows_;
        cols_ = m.cols_;
        allocateMatrixSpace();
    }

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            _matrix[i][j] = m._matrix[i][j];
        }
    }
    return *this;
}

template Matrix<double>& Matrix<double>::operator=(const Matrix<double>& m);

template<typename T>
Matrix<T>& Matrix<T>::operator*=(Matrix<T> const& rhs) {
    if(getCols() == rhs.getRows()){
        Matrix<T> tmp(rows_, rhs.cols_);
        for (int i = 0; i < tmp.rows_; ++i) {
            for (int j = 0; j < tmp.cols_; ++j) {
                for (int k = 0; k < cols_; ++k) {
                    tmp(i,j) += (_matrix[i][k] * rhs.getValue(k,j));
                }
            }
        }
        return (*this = tmp);
    } else {
        throw std::invalid_argument("column dimension of lhs for operation *= is not equal to row dimension of lhs\n");
    }
}

template Matrix<double>& Matrix<double>::operator*=(const Matrix<double>& m);


template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs){
    Matrix<T> result(*this);
    return result *= rhs;
}

template Matrix<double> Matrix<double>::operator*(const Matrix<double>& rhs);

template<typename T>
Matrix<T> Matrix<T>::operator*(const std::vector<T>& rhs){
    if (getCols()==SizeToInt(rhs.size())) {
        // Matrix<T> retvec = Matrix<T>(rows_,1);
        // int m = getRows(), n = getCols();
        // for (int i = 0; i < m; i++) {
        //     for (int k = 0; k < n; k++) {
        //         retvec(i,0) += _matrix[i][k] * rhs[k];
        //     }
        // }
        Matrix<T> retvec = Matrix<T>(rhs);
        return (*this)*retvec;
    } else throw std::invalid_argument("arrayMult vector has not the right size to be multiplied");
}

template Matrix<double> Matrix<double>::operator*(const std::vector<double>& rhs);


template<typename T>
Matrix<T> Matrix<T>::copyAndAddRowsCols(int additionalRows, int additionalCols)const{
    Matrix<T> ret = Matrix<T>(this->getRows()+additionalRows,this->getCols()+additionalCols);
    for (int i = 0 ; i<this->getRows(); i++) {
        for (int j = 0 ; j<this->getRows(); j++) {
            ret(i,j) = getValue(i, j);
        }
    }

    return ret;
}

template Matrix<double> Matrix<double>::copyAndAddRowsCols(int additionalRows, int additionalCols)const;