#include <gtest/gtest.h>
#include <string>
#include <vector>
#include "Computation.h"
#include "Matrix.h"

class ComputationTesting : public ::testing::Test {
    protected:
        void SetUp() override {
            // q0_ remains empty
            c0  = new Computation();
            c1  = new Computation();
            
        }
        void TearDown() override{
            delete c0;
            delete c1;
        }


        std::string thisCellType = "testCell";
        std::vector<double> input= {0.1,0.3,0.5,0.62,0.34,0.87}; 
        Matrix<double> _W=Matrix<double>::createRandom(5, 5); 
        std::vector<std::string> metapathwayNames = {"testCell","testCell2","testCell3","testCell4","testCell5"};

        Computation* c0;       //testing default constructor
        Computation* c1;       //testing general constructor
  
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
