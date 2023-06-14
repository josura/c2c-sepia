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
        }
        void TearDown() override{
            delete c0;
            delete c1;
        }
        ConservationModel* c0;       //testing default constructor
        ConservationModel* c1;       //testing general constructor
  
};

TEST_F(ConservationModelTesting, constructorWorksDefault) {
    EXPECT_EQ(c0->getScaleFunction()(0),0.5);
}

TEST_F(ConservationModelTesting, constructorWorksGeneral) {
    EXPECT_EQ(c1->getScaleFunction()(0),0);
}

// TEST_F(ConservationModelTesting, conservateWorks) {
//     arma::Col<double> input = {1,2,3};
//     arma::Col<double> inputDissipated = {1,2,3};
//     arma::Mat<double> Wstar = {1,2,3};
//     double time = 0;
//     std::vector<double> q = {1,2,3};
//     EXPECT_EQ(c0->conservate(input,inputDissipated,Wstar,time,q)(0),0.5);
// }

