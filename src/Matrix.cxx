#include "Matrix.h"


template<typename T>
void Matrix<T>::allocateMatrixSpace()
{
    _matrix = new T*[rows_];
    for (int i = 0; i < rows_; ++i) {
        _matrix[i] = new T[cols_];
    }
}

template<typename T>
Matrix<T>::Matrix(int rows, int cols) : rows_(rows), cols_(cols)
{
    allocateMatrixSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            _matrix[i][j] = 0;
        }
    }
}

template<typename T>
Matrix<T>::Matrix(T** a, int rows, int cols) : rows_(rows), cols_(cols)
{
    allocateMatrixSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            _matrix[i][j] = a[i][j];
        }
    }
}

template<typename T>
Matrix<T>::Matrix() : rows_(1), cols_(1)
{
    allocateMatrixSpace();
    _matrix[0][0] = 0;
}

template<typename T>
Matrix<T>::~Matrix()
{
    for (int i = 0; i < rows_; ++i) {
        delete[] _matrix[i];
    }
    delete[] _matrix;
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& m) : rows_(m.rows_), cols_(m.cols_)
{
    allocateMatrixSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            _matrix[i][j] = m._matrix[i][j];
        }
    }
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
    if (this == &m) {
        return *this;
    }

    if (rows_ != m.rows_ || cols_ != m.cols_) {
        for (int i = 0; i < rows_; ++i) {
            delete[] _matrix[i];
        }
        delete[] _matrix;

        rows_ = m.rows_;
        cols_ = m.cols_;
        allocateMatrixSpace();
    }

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            _matrix[i][j] = m._matrix[i][j];
        }
    }
    return *this;
}
