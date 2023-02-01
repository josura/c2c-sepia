#include <gtest/gtest.h>
#include <string>
#include <vector>
#include "Computation.h"

class ComputationTesting : public ::testing::Test {
 protected:
  void SetUp() override {
    // q0_ remains empty
    c0  = new Computation();
     
  }
  void TearDown() override{
    delete c0;
  }

  Computation* c0;       //testing default constructor
  
};

TEST_F(ComputationTesting, constructorWorksDefault) {
    EXPECT_EQ(c0->getInput().size(),0);
    EXPECT_EQ(c0->getOutput().size(),0);
    EXPECT_EQ(c0->getLocalCellType(),"");
    EXPECT_EQ(c0->getMetapathway().getNumNodes(),0);
    EXPECT_EQ(c0->getAugmentedMetapathway().getNumNodes(),0);
    EXPECT_EQ(c0->getCellTypes().size(),0);
}

TEST_F(ComputationTesting, constructorWorksGeneral) {
    EXPECT_EQ(c0->getInput().size(),0);
    EXPECT_EQ(c0->getOutput().size(),0);
    EXPECT_EQ(c0->getLocalCellType(),"");
    EXPECT_EQ(c0->getMetapathway().getNumNodes(),0);
    EXPECT_EQ(c0->getAugmentedMetapathway().getNumNodes(),0);
    EXPECT_EQ(c0->getCellTypes().size(),0);
}
