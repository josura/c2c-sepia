#include "Matrix.h"


// helper functions
template<typename T>
void Matrix<T>::allocateMatrixSpace()
{
    _matrix = new T*[rows_];
    for (int i = 0; i < rows_; ++i) {
        _matrix[i] = new T[cols_];
    }
}

template void Matrix<double>::allocateMatrixSpace();


//constructors and destructors
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

template Matrix<double>::Matrix(int rows, int cols);

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

template Matrix<double>::Matrix(double** a, int rows, int cols);

template<typename T>
Matrix<T>::Matrix() : rows_(1), cols_(1)
{
    allocateMatrixSpace();
    _matrix[0][0] = 0;
}

template Matrix<double>::Matrix();

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
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
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

//operators
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

template Matrix<double>& Matrix<double>::operator=(const Matrix<double>& m);

template<typename T>
Matrix<T>& Matrix<T>::operator*=(Matrix<T> const& rhs) {
    if(getCols() == rhs.getRows()){
        for (int r = 0; r < getRows(); r++) {
                  for (int c = 0; c < rhs.getCols(); c++) {
                          //result(r,c) = 0;   //already initialized at 0
                          for (int i = 0; i < rhs.getRows(); i++)
                                  _matrix[r][c] += getValue(r,i) * rhs.getValue(i,c);
                  }
          }
        return *this;
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