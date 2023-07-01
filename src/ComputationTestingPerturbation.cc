#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>
#include "Computation.h"
#include "Matrix.h"
#include "utilities.h"

class ComputationTestingPerturbation : public ::testing::Test {
    protected:
        void SetUp() override {
            c0  = new Computation();
            c1  = new Computation(thisCellType,input,_W,nodesNames);
            
        }
        void TearDown() override{
            //delete c0;
            delete c1;
        }


        std::string thisCellType = "type1";
        std::vector<double> input{1,1,0.1,2}; 
        //Matrix<double> _W=Matrix<double>::createRandom(6, 6); 
        std::vector<double> matrixVector{0,0.2,0.4,0.5,
                                         0.2,0,0.4,0.5,
                                         0.3,0.2,0,0.5,
                                         0.4,0.2,0.4,0
                                         };  // 30 edges initially
        Matrix<double> _W = Matrix<double>(matrixVector,4,4);
        std::vector<std::string> nodesNames{"node1","node2","node3","node4"};

        const std::vector<std::string> types{"type1","type2","type3"};
        const std::vector<std::pair<std::string, std::string>> virtualInputEdges{{"v-in:type2","node2"},
                                                                                 {"v-in:type3","node2"},
                                                                                 {"v-in:type3","node3"},
                                                                                 {"v-in:type3","node1"}
                                                                                 };
        std::vector<double> virtualInputEdgesValues = {0.4,0.5,0.7,0.2};
        std::vector<std::pair<std::string, std::string>> virtualOutputEdges{{"node2","v-out:type2"},
                                                                            {"node4","v-out:type1"}
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
    EXPECT_EQ(c0->gettypes().size(),0);
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
    
    EXPECT_EQ(c1->gettypes().size(),0);
}

TEST_F(ComputationTesting, testingAddingEdges) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types);
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
    EXPECT_EQ(computationTest.gettypes().size(),3);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:type2","node2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:type3","node2"), 0.5);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:type3","node3"), 0.7);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:type3","node1"), 0.2);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("node2","v-out:type2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("node4","v-out:type1"), 0.4);
}

TEST_F(ComputationTesting, testingAddingEdgesUndirected) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types);
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
    EXPECT_EQ(computationTest.gettypes().size(),3);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:type2","node2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:type3","node2"), 0.5);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:type3","node3"), 0.7);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("v-in:type3","node1"), 0.2);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("node2","v-in:type2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("node2","v-in:type3"), 0.5);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("node3","v-in:type3"), 0.7);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("node1","v-in:type3"), 0.2);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("node2","v-out:type2"), 0.4);
    EXPECT_DOUBLE_EQ(augMeta->getEdgeWeight("node4","v-out:type1"), 0.4);
}


TEST_F(ComputationTesting, testingAddingEdgesArmaInitialized) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types);
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
    EXPECT_EQ(computationTest.gettypes().size(),3);
    auto inputArma = computationTest.getInputAugmentedArma();

    std::vector<double> normalizationFactors(computationTest.getAugmentedMetapathway()->getNumNodes(),0);
        for (int i = 0; i < computationTest.getAugmentedMetapathway()->getNumNodes(); i++) {
            for(int j = 0; j < computationTest.getAugmentedMetapathway()->getNumNodes();j++){
                double betaToAdd = std::abs(computationTest.getAugmentedMetapathway()->getEdgeWeight(j,i));
                normalizationFactors[i] += betaToAdd; 
            }
        }
    auto wtransArma = computationTest.getAugmentedMetapathway()->adjMatrix.transpose().normalizeByVectorRow(normalizationFactors).asArmadilloMatrix();
    EXPECT_EQ(inputArma.n_cols, 1);
    EXPECT_EQ(inputArma.n_rows, 12);
    EXPECT_EQ(wtransArma.n_cols, 12);
    EXPECT_EQ(wtransArma.n_rows, 12);
    //edges before
    EXPECT_NEAR(augMeta->getEdgeWeight("v-in:type2","node2"),0.4,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight("v-in:type1","node2"),0,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight("v-in:type3","node2"),0.5,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight(augMeta->getIndexFromName("v-in:type2"),augMeta->getIndexFromName("node2")),0.4,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight(augMeta->getIndexFromName("v-in:type1"),augMeta->getIndexFromName("node2")),0,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight(augMeta->getIndexFromName("v-in:type3"),augMeta->getIndexFromName("node2")),0.5,1e-6);
    //edges after normalization
    EXPECT_NEAR(wtransArma(augMeta->getIndexFromName("node2"),
                         augMeta->getIndexFromName("v-in:type2")),0.210526,1e-6); //0.4 not normalized, inverted since the matrix is transposed
    EXPECT_NEAR(wtransArma(augMeta->getIndexFromName("node2"),
                         augMeta->getIndexFromName("v-in:type3")),0.263157,1e-6); //0.5 not normalized
    EXPECT_NEAR(wtransArma(augMeta->getIndexFromName("node3"),
                         augMeta->getIndexFromName("v-in:type3")),0.259259,1e-6); //0.7 not normalized
    EXPECT_NEAR(wtransArma(augMeta->getIndexFromName("node1"),
                         augMeta->getIndexFromName("v-in:type3")),0.056338,1e-6); //0.2 not normalized
    //control if the change of values is in place (not working) or in copy(working correctly)
    EXPECT_NEAR(augMeta->getEdgeWeight("v-in:type2","node2"),0.4,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight("v-in:type1","node2"),0,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight("v-in:type3","node2"),0.5,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight(augMeta->getIndexFromName("v-in:type2"),augMeta->getIndexFromName("node2")),0.4,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight(augMeta->getIndexFromName("v-in:type1"),augMeta->getIndexFromName("node2")),0,1e-6);
    EXPECT_NEAR(augMeta->getEdgeWeight(augMeta->getIndexFromName("v-in:type3"),augMeta->getIndexFromName("node2")),0.5,1e-6);
}

TEST_F(ComputationTesting, testingAugmentingPathwayNoSelf) {
    Computation computationTest;
    computationTest.assign(*c1);
    
    computationTest.augmentMetapathway(types);
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
    EXPECT_EQ(computationTest.gettypes().size(),3);
}

TEST_F(ComputationTesting, testingAugmentingPathwaySelf) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types,virtualInputEdges,virtualInputEdgesValues,true);
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
    EXPECT_EQ(computationTest.gettypes().size(),4);
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
    computationTest.augmentMetapathway(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
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
    EXPECT_EQ(computationTest.gettypes().size(),3);
}


TEST_F(ComputationTesting, testComputePerturbationSelf){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types,virtualInputEdges,virtualInputEdgesValues,true);
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
    EXPECT_EQ(computationTest.gettypes().size(),4);
}

TEST_F(ComputationTesting, testComputeAugmentedPerturbationSelf){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types,virtualInputEdges,virtualInputEdgesValues,true);
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
    EXPECT_EQ(computationTest.gettypes().size(),4);
}


TEST_F(ComputationTesting, testUpdateInputDefault){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types,virtualInputEdges,virtualInputEdgesValues);
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
    EXPECT_EQ(computationTest.gettypes().size(),3);
}

TEST_F(ComputationTesting, testUpdateInputPerturbation){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types,virtualInputEdges,virtualInputEdgesValues);
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
    EXPECT_EQ(computationTest.gettypes().size(),3);
}


TEST_F(ComputationTesting, testUpdateInputAugmentedDefaultSelf){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types,virtualInputEdges,virtualInputEdgesValues,true);
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
    EXPECT_EQ(computationTest.gettypes().size(),4);
}

TEST_F(ComputationTesting, testUpdateInputAugmentedSelf){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types,virtualInputEdges,virtualInputEdgesValues,true);
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
    EXPECT_EQ(computationTest.gettypes().size(),4);
}


TEST_F(ComputationTesting, testAugmentedVinputZeros){
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathway(types);
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
    EXPECT_EQ(computationTest.gettypes().size(),3);
    ASSERT_TRUE(augMeta->getIndexFromName("v-in:testCell") < 0);
    ASSERT_TRUE(augMeta->getIndexFromName("v-in:type2") >= 0 && augMeta->getIndexFromName("v-in:type2") < 12);
    EXPECT_NEAR(perturbation[augMeta->getIndexFromName("v-in:type2")], 0, 1e-12);
    ASSERT_TRUE(augMeta->getIndexFromName("v-in:type2") >= 0 && augMeta->getIndexFromName("v-in:type1") < 12);
    EXPECT_NEAR(perturbation[augMeta->getIndexFromName("v-in:type1")], 0,1e-12);
    ASSERT_TRUE(augMeta->getIndexFromName("v-in:type2") >= 0 && augMeta->getIndexFromName("v-in:type3") < 12);
    EXPECT_NEAR(perturbation[augMeta->getIndexFromName("v-in:type3")], 0,1e-12);
}

//TESTING IF NODE VALUES FOR v-input nodes are the same as the previous iteration

//TODO TESTING FOR THROWS