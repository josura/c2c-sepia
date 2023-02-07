#include <gtest/gtest.h>
#include <iostream>
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
        std::vector<double> input{0.1,0.3,0.5,0.62,0.34,0.87}; 
        //Matrix<double> _W=Matrix<double>::createRandom(6, 6); 
        std::vector<double> matrixVector{0,0.2,0.4,0.5,0.23,0.67,
                                         0.2,0,0.4,0.5,0.23,0.67,
                                         0.3,0.2,0,0.5,0.23,0.67,
                                         0.43,0.2,0.4,0,0.23,0.67,
                                         0.54,0.2,0.4,0.5,0,0.67,
                                         0.25,0.2,0.4,0.5,0.23,0,};
        Matrix<double> _W = Matrix<double>(matrixVector,6,6);
        std::vector<std::string> metapathwayNames{"testGene1","testGene2","testGene3","testGene4","testGene5","testGene6"};

        const std::vector<std::string> cellTypes{"testCell","testCell2","testCell3","testCell4"};
        const std::vector<std::pair<std::string, std::string>> virtualInputEdges{{"v-in:testCell2","testGene2"},
                                                                                 {"v-in:testCell4","testGene2"},
                                                                                 {"v-in:testCell4","testGene3"},
                                                                                 {"v-in:testCell4","testGene6"}
                                                                                 };
        std::vector<double> virtualInputEdgesValues = {0.4,0.5,0.7,0.2};
        std::vector<std::pair<std::string, std::string>> virtualOutputEdges{{"testGene2","v-out:testCell2"},
                                                                                           {"testGene4","v-out:testCell3"}
                                                                                     };
        std::vector<double> virtualOutputEdgesValues = {0.4,0.4};

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
    EXPECT_EQ(c1->getInputAugmented().size(),0);
    EXPECT_EQ(c1->getOutputAugmented().size(),0);
    EXPECT_EQ(c1->getLocalCellType(),"testCell");
    auto meta = c1->getMetapathway();
    auto augMeta = c1->getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),0);
    
    EXPECT_EQ(c1->getCellTypes().size(),0);
}

TEST_F(ComputationTesting, testingAugmentingPathwayNoSelf) {
    Computation computationTest;
    computationTest.assign(*c1);
    
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),0);
    EXPECT_EQ(computationTest.getInputAugmented().size(),12);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),0);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),15);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),21);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
}

TEST_F(ComputationTesting, testingAugmentingPathwaySelf) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues,true);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),0);
    EXPECT_EQ(computationTest.getInputAugmented().size(),14);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),0);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),15);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),14);
    EXPECT_EQ(augMeta->getNumEdges(),21);
    EXPECT_EQ(computationTest.getCellTypes().size(),4);
}

TEST_F(ComputationTesting, testComputePerturbation){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computePerturbation();
    EXPECT_EQ( perturbation.size(), 6);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),6);
    EXPECT_EQ(computationTest.getInputAugmented().size(),12);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),0);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),15);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),21);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
}

TEST_F(ComputationTesting, testComputeAugmentedPerturbation){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computeAugmentedPerturbation();
    EXPECT_EQ( perturbation.size(), 12);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),0);
    EXPECT_EQ(computationTest.getInputAugmented().size(),12);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),12);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),15);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),21);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
}
