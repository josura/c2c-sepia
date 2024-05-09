#include "gtest/gtest.h"
#include <iostream>
#include <string>
#include <vector>
#include "ConservationModel.h"


class ConservationModelTesting : public ::testing::Test {
    protected:
        void SetUp() override {
            c0 = new ConservationModel();
            c1 = new ConservationModel([](double time)->double{return 0;});

            Wstar_oneEdge(0,1) = 1;

            Wstar_threeEdges(0,1) = 1;
            Wstar_threeEdges(1,2) = 1;
            Wstar_threeEdges(2,0) = 1;
        }
        void TearDown() override{
            delete c0;
            delete c1;
        }
        ConservationModel* c0;       //testing default constructor
        ConservationModel* c1;       //testing general constructor
        arma::Col<double> input = {1,2,3};
        arma::Col<double> inputDissipated = {1,2,3};
        arma::dmat Wstar_zeros = arma::zeros<arma::dmat>(3,3);
        arma::dmat Wstar_oneEdge = arma::zeros<arma::dmat>(3,3);
        arma::dmat Wstar_threeEdges = arma::zeros<arma::dmat>(3,3);
        std::vector<double> q = {1,1,1};
  
};

TEST_F(ConservationModelTesting, constructorWorksDefault) {
    EXPECT_DOUBLE_EQ(c0->getScaleFunction()(0),0.5);
}

TEST_F(ConservationModelTesting, constructorWorksGeneral) {
    EXPECT_DOUBLE_EQ(c1->getScaleFunction()(0),0);
}

TEST_F(ConservationModelTesting, conservateWorksZeroEdges) {
    double time = 0;
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_zeros,time,q)(0),1);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_zeros,time,q)(1),2);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_zeros,time,q)(2),3);
}

TEST_F(ConservationModelTesting, conservateWorksZeroScaleFunctionZeroEdges) {
    double time = 0;
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_zeros,time)(0),1.0);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_zeros,time)(1),2.0);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_zeros,time)(2),3.0);
    EXPECT_DOUBLE_EQ(c1->conservate(input,inputDissipated,Wstar_zeros,time,q)(0),1);
    EXPECT_DOUBLE_EQ(c1->conservate(input,inputDissipated,Wstar_zeros,time,q)(1),2);
    EXPECT_DOUBLE_EQ(c1->conservate(input,inputDissipated,Wstar_zeros,time,q)(2),3);
}

TEST_F(ConservationModelTesting, conservateWorksOneEdge) {
    double time = 0;

    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_oneEdge,time)(0),0.5);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_oneEdge,time)(1),2.0);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_oneEdge,time)(2),3.0);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_oneEdge,time,q)(0),0.5);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_oneEdge,time,q)(1),2.0);
    EXPECT_DOUBLE_EQ(c0->conservate(input,inputDissipated,Wstar_oneEdge,time,q)(2),3.0);
}

