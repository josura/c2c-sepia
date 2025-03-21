#include "Matrix.h"
#include "utilities.h"
#include <iostream>
#include <stdexcept>
#include <vector>


// helper functions
template<typename T>
void Matrix<T>::allocateMatrixSpace()
{
    if(_matrix) {
        delete[] _matrix;
        _matrix=nullptr;
    }
    int totalLength = rows_ * cols_;
    _matrix = new T[totalLength]; //columns first
    
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            _matrix[i*cols_ + j] = 0;
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
            _matrix[i*cols_ + j] = a[i][j];
        }
    }
}

template Matrix<double>::Matrix(double** a, int rows, int cols);

template<typename T>
Matrix<T>::Matrix() : rows_(1), cols_(1)
{
    allocateMatrixSpace();
    _matrix[0] = 0;
}

template Matrix<double>::Matrix();

template<typename T>
Matrix<T>::Matrix(const std::vector<T>& vec, uint nrows, uint ncols):rows_(vec.size()),cols_(ncols){
    bool isAvector = false;
    if((nrows == vec.size() || nrows == 0) && ncols == 1){
        isAvector=true;
    }
    else if( (nrows != vec.size() && (approximatelyEqual(vec.size()/(double)nrows, (double)ncols, std::numeric_limits<T>::epsilon())))){
        rows_ = nrows;
    } else {
        std::cerr << "[ERROR] Matrix<T>::Matrix(vec,nrows=0,ncols=1): the number of resulting columns from the division of rows is not equal to the one passed in the ncols parameter: vec.size="<< vec.size() << " nrows=" << nrows << " ncols="<<ncols; 
        throw std::invalid_argument("[ERROR] Matrix<T>::Matrix(vec,nrows=0,ncols=1): the number of resulting columns from the division of rows is not equal to the one passed in the ncols parameter");
    }
    allocateMatrixSpace();
    if(isAvector){
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_;j++) {
                int index = i * cols_ + j; 
                _matrix[index] = vec[index];
            }
        }
    }else {
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < cols_;j++) {
                int index = i * cols_ + j; 
                _matrix[index] = vec[index];
            }
        }
    }
    
}

template Matrix<double>::Matrix(const std::vector<double>& vec,uint nrows, uint ncols);


template<typename T>
Matrix<T>::~Matrix()
{
    if(_matrix) {
        delete[] _matrix;
        _matrix=nullptr;
    }
}

template Matrix<double>::~Matrix();

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& m) : rows_(m.rows_), cols_(m.cols_)
{
    allocateMatrixSpace();
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            _matrix[i*cols_ + j] = m.getValue(i,j);
        }
    }
}

template Matrix<double>::Matrix(const Matrix<double>& m);

//static methods
// generate random matrix, values ranging from -DOUBLE_MAX to DOUBLE_MAX
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
        if(_matrix){
            delete[] _matrix;
            _matrix = nullptr;
        }

        rows_ = m.rows_;
        cols_ = m.cols_;
        allocateMatrixSpace();
    }

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            _matrix[i*cols_ + j] = m.getValue(i,j);
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
                    tmp(i,j) += (getValue(i, k) * rhs.getValue(k,j));
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
Matrix<T>& Matrix<T>::operator*=(T rhs) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            _matrix[i * cols_ + j] *= rhs;
        }
    }
    return *this;
}

template Matrix<double>& Matrix<double>::operator*=(double m);


template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m)
{
    if(getCols() == m.getRows()){
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                _matrix[i * cols_ + j] += m.getValue(i,j);
            }
        }
        return *this;
    } else {
        throw std::invalid_argument("column dimension of lhs for operation += is not equal to row dimension of lhs\n");
    }
}

template Matrix<double>& Matrix<double>::operator+=(const Matrix<double>& m);


template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m)
{
    if(getCols() == m.getRows()){
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                _matrix[i * cols_ + j] -= m.getValue(i,j);
            }
        }
        return *this;
    } else {
        throw std::invalid_argument("column dimension of lhs for operation -= is not equal to row dimension of lhs\n");
    }
}

template Matrix<double>& Matrix<double>::operator-=(const Matrix<double>& m);

template<typename T>
Matrix<T>& Matrix<T>::operator/=(T num)
{
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            _matrix[i * cols_ + j] /= num;
        }
    }
    return *this;
}

template Matrix<double>& Matrix<double>::operator/=(double m);


template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs)const{
    Matrix<T> result(*this);
    return result *= rhs;
}

template Matrix<double> Matrix<double>::operator*(const Matrix<double>& rhs)const;

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs)const{
    Matrix<T> result(*this);
    return result -= rhs;
}

template Matrix<double> Matrix<double>::operator-(const Matrix<double>& rhs)const;

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& rhs)const{
    Matrix<T> result(*this);
    return result += rhs;
}

template Matrix<double> Matrix<double>::operator+(const Matrix<double>& rhs)const;


template<typename T>
Matrix<T> Matrix<T>::operator*(const std::vector<T>& rhs)const{
    if (getCols()==SizeToInt(rhs.size())) {
        // Matrix<T> retvec = Matrix<T>(rows_,1);
        // int m = getRows(), n = getCols();
        // for (int i = 0; i < m; i++) {
        //     for (int k = 0; k < n; k++) {
        //         retvec(i,0) += getValue(i, k) * rhs[k];
        //     }
        // }
        Matrix<T> retvec = Matrix<T>(rhs);
        return (*this)*retvec;
    } else throw std::invalid_argument("arrayMult vector has not the right size to be multiplied");
}

template Matrix<double> Matrix<double>::operator*(const std::vector<double>& rhs)const;



template<typename T>
Matrix<T> Matrix<T>::transpose()const
{
    Matrix ret(cols_, rows_);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            ret(j,i) = getValue(i, j);
        }
    }
    return ret;
}

template Matrix<double> Matrix<double>::transpose()const;


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

template<typename T>
std::vector<T> Matrix<T>::asVector()const{
    if(this->isVector()){
        std::vector<T> ret(rows_,0);
        for (int i = 0; i<rows_; i++) {
            ret[i]=getValue(i,0);
        }
        return ret;
    } else {
        std::cerr << "[ERROR] Matrix::asVector: the matrix is not a vector (1 column, n rows)";
        std::cerr << "rows_=" << rows_ << " cols_=" << cols_ << std::endl;
        throw std::domain_error("[ERROR] Matrix::asVector: the matrix is not a vector (1 column, n rows)");
    }
}

template std::vector<double> Matrix<double>::asVector()const;


template<typename T>
Matrix<T> Matrix<T>::concatenateRight(const Matrix<T>& rhs)const{
    if(rhs.getRows()==rows_){
        Matrix<T> ret = this->copyAndAddRowsCols(0, rhs.cols_);
        for (int i = 0; i<rows_; i++) {
            for (int j = 0; j<rhs.getCols(); j++) {
                ret(i,j+cols_) = rhs.getValue(i,j);    
            }
            
        }
        return ret;
    } else {
        throw std::invalid_argument("[ERROR] Matrix::concatenateRight: rhs is not of the same size(rows)");
    }
}


template Matrix<double> Matrix<double>::concatenateRight(const Matrix<double>& rhs)const;


//minor
template<typename T>
Matrix<T> Matrix<T>::getMinor(const Matrix<T>& A,int p, int q,int n) {
    int i = 0, j = 0;
    Matrix minor(n-1,n-1);
    for (int row = 0; row < n; row++){
        for (int col = 0; col<n; col++){
            if (row != p && col != q){
                minor(i,j++) = A.getValue(row,col);
                if (j == n - 1){
                    j = 0;
                    i++;
                }
            }
        }
    }
    return minor;
}

template Matrix<double> Matrix<double>::getMinor(const Matrix<double>& A,int p, int q,int n);

//determinant
template<typename T>
T Matrix<T>::determinant(const Matrix<T>& A){
    T deter=0;
    int n = A.getRows();
    switch (n) {
        case (1):
            return A.getValue(0, 0);
        case (2):
            return((A.getValue(0,0)*A.getValue(1,1))-(A.getValue(0,1)*A.getValue(1,0))); 
        case (3):
            return((A.getValue(0,0)*A.getValue(1,1)*A.getValue(2,2))+(A.getValue(1,0)*A.getValue(1,2)*A.getValue(2,0))+(A.getValue(1,0)*A.getValue(2,1)*A.getValue(0,2))-(A.getValue(0,2)*A.getValue(1,1)*A.getValue(2,0))-(A.getValue(0,1)*A.getValue(1,0)*A.getValue(2,2))-(A.getValue(1,2)*A.getValue(2,1)*A.getValue(0,0)));
    }
    Matrix tempCofactor(n,0.0); 
    int sign = 1;  // To store sign multiplier 
     // first row fixed
    for (int f = 0; f < n; f++){
        tempCofactor=getMinor(A,0, f, n);
        deter += sign * A.getValue(0,f) * determinant(tempCofactor);         //cofactor =sign * A[0][f] * determinant(tempCofactor, n - 1)
        sign = -sign;
    } 
    return deter;
}

template double Matrix<double>::determinant(const Matrix<double>& A);

template<typename T>
T Matrix<T>::determinant()const{
    T deter=0;
    switch (rows_) {
        case (1):
            return getValue(0, 0);
        case (2):
            return((getValue(0, 0)*getValue(1, 1))-(getValue(0, 1)*getValue(1, 0))); 
        case (3):
            return((getValue(0, 0)*getValue(1, 1)*getValue(2, 2))+(getValue(1, 0)*getValue(1, 2)*getValue(2, 0))+(getValue(1, 0)*getValue(2, 1)*getValue(0, 2))-(getValue(0, 2)*getValue(1, 1)*getValue(2, 0))-(getValue(0, 1)*getValue(1, 0)*getValue(2, 2))-(getValue(1, 2)*getValue(2, 1)*getValue(0, 0)));
    }
    Matrix tempCofactor(rows_,rows_); 
    int sign = 1;  // To store sign multiplier 
     // first row fixed
    for (int f = 0; f < rows_; f++){
    	Matrix minor=*this;
        tempCofactor=getMinor(minor, 0, f, rows_);
        deter += sign * getValue(0, f) * determinant(tempCofactor);         
        sign = -sign;
    } 
    return deter;
}

template double Matrix<double>::determinant()const;

template<typename T>
Matrix<T> Matrix<T>::inverse(){
    auto ret = Matrix<T>();
    //TODO choose between Gaussian elimination method and adjugate method, also controls on determinant after reduction to triangular form?
    return ret;
}

template Matrix<double> Matrix<double>::inverse();

template<typename T>
arma::Mat<T> Matrix<T>::asArmadilloMatrix()const{
    arma::Mat<T> ret = arma::Mat<T>(getRows(),getCols());  //filled with zeros
    for (int i =0 ; i<getRows(); i++) {
        for (int j =0 ; j<getCols(); j++) {
            ret(i,j) = getValue(i, j);
        }
    }
    return ret;
    //return arma::Mat<T>(_matrix,getRows(),getCols());
}

template arma::Mat<double> Matrix<double>::asArmadilloMatrix()const;


template<typename T>
arma::Col<T> Matrix<T>::asArmadilloColumnVector()const{
    return arma::Col<T>(asVector());
}

template arma::Col<double> Matrix<double>::asArmadilloColumnVector()const;


template<typename T>
arma::Row<T> Matrix<T>::asArmadilloRowVector()const{
    return arma::Row<T>(asVector());
}

template arma::Row<double> Matrix<double>::asArmadilloRowVector()const;


template<typename T>
Matrix<T>& Matrix<T>::normalizeByVectorColumn(const std::vector<double>& normVector){    
    if (SizeToInt( normVector.size())< cols_) {
        throw std::invalid_argument("[ERROR] Matrix<T>::normalizeByVectorColumn: normalization vector size is less than the number of columns");
    }
    for (int i =0 ; i<getRows(); i++) {
        for (int j =0 ; j<getCols(); j++) {
            _matrix[i * cols_ + j] = getValue(i, j) / (normVector[j] + 1e-20);
        }
    }
    return *this;
}

template Matrix<double>& Matrix<double>::normalizeByVectorColumn(const std::vector<double> &normVector);

template<typename T>
Matrix<T>& Matrix<T>::normalizeByVectorRow(const std::vector<double>& normVector){
    if (SizeToInt( normVector.size())< rows_) {
        throw std::invalid_argument("[ERROR] Matrix<T>::normalizeByVectorRow: normalization vector size is less than the number of rows");
    }
    for (int i =0 ; i<getRows(); i++) {
        for (int j =0 ; j<getCols(); j++) {
            _matrix[i * cols_ + j] = getValue(i, j) / (normVector[i] + 1e-20) ;
        }
    }
    return *this;
}

template Matrix<double>& Matrix<double>::normalizeByVectorRow(const std::vector<double> &normVector);

template<typename T>
void Matrix<T>::printMatrix()const{
    std::cout << "Matrix " << rows_ << "x" << cols_ << std::endl;
    for (int i =0 ; i<getRows(); i++) {
        for (int j =0 ; j<getCols(); j++) {
            std::cout << getValue(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

template void Matrix<double>::printMatrix()const;