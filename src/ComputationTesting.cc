#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>
#include "Computation.h"
#include "Matrix.h"

class ComputationTesting : public ::testing::Test {
    protected:
        void SetUp() override {
            c0  = new Computation();
            c1  = new Computation(thisCellType,input,_W,metapathwayNames);
            
        }
        void TearDown() override{
            //delete c0;
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
                                         0.25,0.2,0.4,0.5,0.23,0,};  // 30 edges initially
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),0);
    
    EXPECT_EQ(c1->getCellTypes().size(),0);
}

TEST_F(ComputationTesting, testingAddingEdges) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:testCell2","testGene2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:testCell4","testGene2"), 0.5);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:testCell4","testGene3"), 0.7);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:testCell4","testGene6"), 0.2);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("testGene2","v-out:testCell2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("testGene4","v-out:testCell3"), 0.4);
}

TEST_F(ComputationTesting, testingAddingEdgesUndirected) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues,true);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues,true);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),0);
    EXPECT_EQ(computationTest.getInputAugmented().size(),12);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),0);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),42);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:testCell2","testGene2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:testCell4","testGene2"), 0.5);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:testCell4","testGene3"), 0.7);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:testCell4","testGene6"), 0.2);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("testGene2","v-in:testCell2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("testGene2","v-in:testCell4"), 0.5);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("testGene3","v-in:testCell4"), 0.7);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("testGene6","v-in:testCell4"), 0.2);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("testGene2","v-out:testCell2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("testGene4","v-out:testCell3"), 0.4);
}


TEST_F(ComputationTesting, testingAddingEdgesArmaInitialized) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
    auto inputArma = computationTest.getInputAugmentedArma();
    auto identityArma = computationTest.getIdentityAugmentedArma();
    auto wtransArma = computationTest.getWtransAugmentedArma();
    EXPECT_EQ(inputArma.n_cols, 1);
    EXPECT_EQ(inputArma.n_rows, 12);
    EXPECT_EQ(identityArma.n_rows, 12);
    EXPECT_EQ(identityArma.n_cols, 12);
    EXPECT_EQ(wtransArma.n_cols, 12);
    EXPECT_EQ(wtransArma.n_rows, 12);
    EXPECT_EQ(wtransArma(augMeta->getIndexFromName("testGene2"),
                         augMeta->getIndexFromName("v-in:testCell2")),0.4); //inverted since the matrix is transposed
    EXPECT_EQ(wtransArma(augMeta->getIndexFromName("testGene2"),
                         augMeta->getIndexFromName("v-in:testCell4")),0.5);
    EXPECT_EQ(wtransArma(augMeta->getIndexFromName("testGene3"),
                         augMeta->getIndexFromName("v-in:testCell4")),0.7);
    EXPECT_EQ(wtransArma(augMeta->getIndexFromName("testGene6"),
                         augMeta->getIndexFromName("v-in:testCell4")),0.2);
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),36);
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),14);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),4);
}

TEST_F(ComputationTesting, testComputePerturbation){
    Computation computationTest;
    computationTest.assign(*c1);
    auto perturbation = computationTest.computePerturbation();
    EXPECT_EQ( perturbation.size(), 6);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),6);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),30);
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
}


TEST_F(ComputationTesting, testComputePerturbationSelf){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues,true);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computePerturbation();
    EXPECT_EQ( perturbation.size(), 6);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),6);
    EXPECT_EQ(computationTest.getInputAugmented().size(),14);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),0);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),14);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),4);
}

TEST_F(ComputationTesting, testComputeAugmentedPerturbationSelf){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues,true);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computeAugmentedPerturbation();
    EXPECT_EQ( perturbation.size(), 14);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),0);
    EXPECT_EQ(computationTest.getInputAugmented().size(),14);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),14);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),14);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),4);
}


TEST_F(ComputationTesting, testUpdateInputDefault){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computePerturbation();
    computationTest.updateInput();
    auto computationInputAfter = computationTest.getInput();
    EXPECT_EQ(computationInputAfter.size(), perturbation.size());
    for (uint i = 0; i< computationInputAfter.size(); i++) {
        EXPECT_DOUBLE_EQ(perturbation[i], computationInputAfter[i]);
    }
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
}

TEST_F(ComputationTesting, testUpdateInputPerturbation){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computePerturbation();
    computationTest.updateInput(perturbation);
    auto computationInputAfter = computationTest.getInput();
    EXPECT_EQ(computationInputAfter.size(), perturbation.size());
    for (uint i = 0; i< computationInputAfter.size(); i++) {
        EXPECT_DOUBLE_EQ(perturbation[i], computationInputAfter[i]);
    }
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
}


TEST_F(ComputationTesting, testUpdateInputAugmentedDefaultSelf){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues,true);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computeAugmentedPerturbation();
    auto nothing = std::vector<double>();
    computationTest.updateInput(nothing,true);
    auto computationInputAfter = computationTest.getInputAugmented();
    EXPECT_EQ(computationInputAfter.size(), perturbation.size());
    for (uint i = 0; i< computationInputAfter.size(); i++) {
        EXPECT_DOUBLE_EQ(perturbation[i], computationInputAfter[i]);
    }
    EXPECT_EQ( perturbation.size(), 14);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),0);
    EXPECT_EQ(computationTest.getInputAugmented().size(),14);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),14);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),14);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),4);
}

TEST_F(ComputationTesting, testUpdateInputAugmentedSelf){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes,virtualInputEdges,virtualInputEdgesValues,true);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computeAugmentedPerturbation();
    computationTest.updateInput(perturbation,true);
    auto computationInputAfter = computationTest.getInputAugmented();
    EXPECT_EQ(computationInputAfter.size(), perturbation.size());
    for (uint i = 0; i< computationInputAfter.size(); i++) {
        EXPECT_DOUBLE_EQ(perturbation[i], computationInputAfter[i]);
    }
    EXPECT_EQ( perturbation.size(), 14);
    EXPECT_EQ(computationTest.getInput().size(),6);
    EXPECT_EQ(computationTest.getOutput().size(),0);
    EXPECT_EQ(computationTest.getInputAugmented().size(),14);
    EXPECT_EQ(computationTest.getOutputAugmented().size(),14);
    EXPECT_EQ(computationTest.getLocalCellType(),"testCell");
    auto meta = computationTest.getMetapathway();
    auto augMeta = computationTest.getAugmentedMetapathway();
    ASSERT_TRUE(meta != nullptr);
    EXPECT_EQ(meta->getNumNodes(),6);
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),14);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),4);
}


TEST_F(ComputationTesting, testAugmentedVinputZeros){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(cellTypes);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = computationTest.computeAugmentedPerturbation();
    computationTest.updateInput(perturbation,true);
    auto computationInputAfter = computationTest.getInputAugmented();
    EXPECT_EQ(computationInputAfter.size(), perturbation.size());
    for (uint i = 0; i< computationInputAfter.size(); i++) {
        EXPECT_DOUBLE_EQ(perturbation[i], computationInputAfter[i]);
    }
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
    EXPECT_EQ(meta->getNumEdges(),30);
    ASSERT_TRUE(augMeta != nullptr);
    EXPECT_EQ(augMeta->getNumNodes(),12);
    EXPECT_EQ(augMeta->getNumEdges(),36);
    EXPECT_EQ(computationTest.getCellTypes().size(),3);
    ASSERT_TRUE(augMeta->getIndexFromName("v-in:testCell") < 0);
    ASSERT_TRUE(augMeta->getIndexFromName("v-in:testCell2") >= 0 && augMeta->getIndexFromName("v-in:testCell2") < 12);
    EXPECT_NEAR(perturbation[augMeta->getIndexFromName("v-in:testCell2")], 0, 1e-12);
    ASSERT_TRUE(augMeta->getIndexFromName("v-in:testCell2") >= 0 && augMeta->getIndexFromName("v-in:testCell3") < 12);
    EXPECT_NEAR(perturbation[augMeta->getIndexFromName("v-in:testCell3")], 0,1e-12);
    ASSERT_TRUE(augMeta->getIndexFromName("v-in:testCell2") >= 0 && augMeta->getIndexFromName("v-in:testCell4") < 12);
    EXPECT_NEAR(perturbation[augMeta->getIndexFromName("v-in:testCell4")], 0,1e-12);
}

//TESTING IF NODE VALUES FOR v-input nodes are the same as the previous iteration

//TODO TESTING FOR THROWS