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
            c1  = new Computation(thisCellType,input,_W,metapathwayNames);
            
        }
        void TearDown() override{
            delete c0;
            delete c1;
        }


        std::string thisCellType = "testCell";
        std::vector<double> input= {0.1,0.3,0.5,0.62,0.34,0.87}; 
        Matrix<double> _W=Matrix<double>::createRandom(6, 6); 
        std::vector<std::string> metapathwayNames = {"testGene1","testGene2","testGene3","testGene4","testGene5","testGene6"};

        Computation* c0;       //testing default constructor
        Computation* c1;       //testing general constructor
  
};

TEST_F(ComputationTesting, constructorWorksDefault) {
    EXPECT_EQ(c0->getInput().size(),0);
    EXPECT_EQ(c0->getOutput().size(),0);
    EXPECT_EQ(c0->getLocalCellType(),"");
    auto meta = c0->getMetapathway();
    auto augMeta = c0->getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),0);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),0);
    EXPECT_EQ(c0->getCellTypes().size(),0);
}

TEST_F(ComputationTesting, constructorWorksGeneral) {
    EXPECT_EQ(c1->getInput().size(),6);
    EXPECT_EQ(c1->getOutput().size(),0);
    EXPECT_EQ(c1->getLocalCellType(),"testCell");
    auto meta = c1->getMetapathway();
    auto augMeta = c1->getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),0);
    
    EXPECT_EQ(c1->getCellTypes().size(),0);
}
