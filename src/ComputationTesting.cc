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
    
}
