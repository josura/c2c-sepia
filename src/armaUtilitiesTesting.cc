#include <gtest/gtest.h>
#include <armadillo>
#include <iostream>
#include <string>
#include <vector>
#include "armaUtilities.h"

class armaUtilitiesTesting : public ::testing::Test {
    protected:
        void SetUp() override {
            //initialize resources
            matr(0,0) = 1;
            matr(0,1) = 2;
            matr(0,2) = 3;
            matr(1,0) = 4;
            matr(1,1) = 5;
            matr(1,2) = 6;
            
        }
        void TearDown() override{
            //delete resources
        }
        // protected variables here
        arma::Mat<double> matr = arma::Mat<double>(2,3);
        arma::Mat<double> zeroMatr = arma::zeros<arma::Mat<double>>(2,3);
        arma::Mat<double> identity = arma::eye<arma::Mat<double>>(3,3);

  
};

TEST_F(armaUtilitiesTesting, matrixInitializationWorks) {
    EXPECT_EQ(matr(0,0), 1);
    EXPECT_EQ(matr(0,1), 2);
    EXPECT_EQ(matr(0,2), 3);
    EXPECT_EQ(matr(1,0), 4);
    EXPECT_EQ(matr(1,1), 5);
    EXPECT_EQ(matr(1,2), 6);

    EXPECT_EQ(zeroMatr(0,0), 0);
    EXPECT_EQ(zeroMatr(0,1), 0);
    EXPECT_EQ(zeroMatr(0,2), 0);
    EXPECT_EQ(zeroMatr(1,0), 0);
    EXPECT_EQ(zeroMatr(1,1), 0);
    EXPECT_EQ(zeroMatr(1,2), 0);

    EXPECT_EQ(identity(0,0), 1);
    EXPECT_EQ(identity(0,1), 0);
    EXPECT_EQ(identity(0,2), 0);
    EXPECT_EQ(identity(1,0), 0);
    EXPECT_EQ(identity(1,1), 1);
    EXPECT_EQ(identity(1,2), 0);
    EXPECT_EQ(identity(2,0), 0);
    EXPECT_EQ(identity(2,1), 0);
    EXPECT_EQ(identity(2,2), 1);
}