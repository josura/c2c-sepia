#include <gtest/gtest.h>
#include "Matrix.h"

class MatrixTesting : public ::testing::Test {
 protected:
  void SetUp() override {
    // q0_ remains empty
    m0_  = new Matrix<double>();
    m1_ = new Matrix<double>(5, 10);
    double** tmp = new double*[10];
    for(int i = 0;i < 10; i++){
      tmp[i] = new double[12];
    }
    m2_ = new Matrix<double>(tmp, 10, 12);
    m3_ = new Matrix<double>(*m2_);

    
    //*g2_ = *g1_;   //Problem with the equal operator(assignment)
     
  }
  void TearDown() override{
    delete m0_;
    delete m1_;
    delete m2_;
    delete m3_;
  }

  // void TearDown() override {}

  Matrix<double>* m0_;       //testing default constructor
  Matrix<double>* m1_;    //testing rows and cols constructor
  Matrix<double>* m2_;  //testing passing another matrix as array input to constructor
  Matrix<double>* m3_;  //testing passing another matrix as input to constructor (= operator)
};

TEST_F(MatrixTesting, constructorWorksDefault) {
  EXPECT_EQ(m0_->getCols(), 1);
  EXPECT_EQ(m0_->getRows(), 1);
  EXPECT_FLOAT_EQ(m0_->getValue(0,0), 0);
}


TEST_F(MatrixTesting, constructorWorksRowsAndCols) {
  EXPECT_EQ(m1_->getCols(), 10);
  EXPECT_EQ(m1_->getRows(), 5);
  for (int i = 0; i<m1_->getRows(); i++) {
    for (int j = 0; j<m1_->getCols(); j++) {
      EXPECT_FLOAT_EQ(m1_->getValue(i,j), 0);
    }
  }
}

TEST_F(MatrixTesting, constructorWorksPassingArray) {
  EXPECT_EQ(m2_->getCols(), 12);
  EXPECT_EQ(m2_->getRows(), 10);
  for (int i = 0; i<m2_->getRows(); i++) {
    for (int j = 0; j<m2_->getCols(); j++) {
      EXPECT_FLOAT_EQ(m2_->getValue(i,j), 0);
    }
  }
}

TEST_F(MatrixTesting, constructorWorksPassingMatrix) {
  EXPECT_EQ(m3_->getCols(), 12);
  EXPECT_EQ(m3_->getRows(), 10);
  for (int i = 0; i<m3_->getRows(); i++) {
    for (int j = 0; j<m3_->getCols(); j++) {
      EXPECT_FLOAT_EQ(m3_->getValue(i,j), 0);
    }
  }
}

TEST_F(MatrixTesting,multiplicationControlDimensions){
  Matrix<double> matrixres = (*m1_) * (*m2_);
  EXPECT_EQ(matrixres.getRows(), 5);
  EXPECT_EQ(matrixres.getCols(), 12);
}

TEST_F(MatrixTesting,multiplicationControlResults){
  Matrix<double> matrixres = (*m1_) * (*m2_);
  for (int i=0;i<matrixres.getRows() ;i++) {
      for (int j=0;j<matrixres.getCols() ;j++) {
          EXPECT_FLOAT_EQ( matrixres.getValue(i,j),0);
      }
  }
}

/*TEST_F(GraphTesting, assignmentWorks) {
  EXPECT_EQ(g2_->getNumNodes(), 4);
  EXPECT_EQ(g2_->getNumEdges(), 2);
}*/
