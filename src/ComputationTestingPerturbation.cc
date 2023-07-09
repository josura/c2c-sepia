#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>
#include "Computation.h"
#include "ConservationModel.h"
#include "DissipationModel.h"
#include "DissipationModelScaled.h"
#include "Matrix.h"
#include "utilities.h"

class ComputationTestingPerturbation : public ::testing::Test {
    protected:
        void SetUp() override {
            c1  = new Computation(thisType,input,_W,nodesNames);
            
        }
        void TearDown() override{
            //delete c0;
            delete c1;
        }


        std::string thisType = "type1";
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
        std::vector<double> virtualInputEdgesValues = {1,1,1,1};
        std::vector<std::pair<std::string, std::string>> virtualOutputEdges{{"node2","v-out:type2"},
                                                                            {"node4","v-out:type3"}
                                                                                     };
        std::vector<double> virtualOutputEdgesValues = {1,1};

        Computation* c1;       //testing general constructor

        DissipationModel* dms = new DissipationModelScaled( [](double x){return 0;} );
        DissipationModel* dms2 = new DissipationModelScaled(); //default constructor uses the functions that returns 0.5 for every iteration

        ConservationModel* cms = new ConservationModel( [](double x){return 0;} );
        ConservationModel* cms2 = new ConservationModel(); //default constructor uses the functions that returns 0.5 for every iteration
  
};

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationNoConservationDefaultNoSaturation) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms);
    computationTest.setConservationModel(cms);
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,false);
    std::vector<double> expected{3.256667,2.736364,2.924242,5.593939,0,0,1.303030,2.796970};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-4);
    }
}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationNoConservationDefaultSaturationSingle) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms);
    computationTest.setConservationModel(cms);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,true,saturationVector);
    std::vector<double> expected{0.5,0.5,0.5,0.5,0,0,0.5,0.5};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-2);
    }
}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectDissipationDefaultNoConservationNoSaturation) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms2);
    computationTest.setConservationModel(cms);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,false);
    std::vector<double> expected{3.256667/2,2.736364/2,2.924242/2,5.593939/2,0,0,1.303030/2,2.796970/2};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-2);
    }
}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectDissipationDefaultNoConservationSaturationsingle) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms2);
    computationTest.setConservationModel(cms);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,true, saturationVector);
    std::vector<double> expected{0.5,0.5,0.5,0.5,0,0,0.5,0.5};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-1);
    }
}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationConservationDefaultQemptyNoSaturation) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms);
    computationTest.setConservationModel(cms2);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,false);
    std::vector<double> expected{3.256667-(input[0]/2),2.736364-(input[1]/2),2.924242-(input[2]/2),5.593939-(input[3]/2),0,0,1.303030,2.796970};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-2);
    }
}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationConservationDefaultQemptySaturationsingle) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms);
    computationTest.setConservationModel(cms2);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,true, saturationVector);
    //std::vector<double> expected{3.256667-(input[0]/2),2.736364-(input[1]/2),2.924242-(input[2]/2),5.593939-(input[3]/2),0,0,1.303030,2.796970};
    std::vector<double> expected{0.5,0.5,0.5,0.5,0,0,0.5,0.5};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-2);
    }
}


TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationConservationQpassedNoSaturation) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms);
    computationTest.setConservationModel(cms2);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    std::vector<double> qVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,false,std::vector<double>(),qVector);
    // std::vector<double> expected{3.256667-(input[0]/2),2.736364-(input[1]/2),2.924242-(input[2]/2),5.593939-(input[3]/2),0,0,1.303030,2.796970};
    std::vector<double> expected{3.256667-(input[0]/4),2.736364-(input[1]/4),2.924242-(input[2]/4),5.593939-(input[3]/4),0,0,1.303030,2.796970};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-2);
    }
}


TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationConservationQpassedSaturationdefault) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms);
    computationTest.setConservationModel(cms2);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    std::vector<double> qVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,true,std::vector<double>(),qVector);
    // std::vector<double> expected{3.256667-(input[0]/2),2.736364-(input[1]/2),2.924242-(input[2]/2),5.593939-(input[3]/2),0,0,1.303030,2.796970};
    // std::vector<double> expected{3.256667-(input[0]/4),2.736364-(input[1]/4),2.924242-(input[2]/4),5.593939-(input[3]/4),0,0,1.303030,2.796970};
    std::vector<double> expected{1,1,1,1,0,0,1,1};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],2e-1);
    }
}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationConservationQpassedSaturationsingle) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms);
    computationTest.setConservationModel(cms2);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    std::vector<double> qVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,true,saturationVector,qVector);
    // std::vector<double> expected{3.256667-(input[0]/2),2.736364-(input[1]/2),2.924242-(input[2]/2),5.593939-(input[3]/2),0,0,1.303030,2.796970};
    // std::vector<double> expected{3.256667-(input[0]/4),2.736364-(input[1]/4),2.924242-(input[2]/4),5.593939-(input[3]/4),0,0,1.303030,2.796970};
    std::vector<double> expected{0.5,0.5,0.5,0.5,0,0,0.5,0.5};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-2);
    }
}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectDissipationConservationQemptyNoSaturation) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms2);
    computationTest.setConservationModel(cms2);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    std::vector<double> qVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,false);
    std::vector<double> expected{3.256667/2-(input[0]/4),2.736364/2-(input[1]/4),2.924242/2-(input[2]/4),5.593939/2-(input[3]/4),0,0,1.303030/2,2.796970/2};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-2);
    }
}


TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectDissipationConservationQpassedNoSaturation) {
    Computation computationTest;
    computationTest.assign(*c1);
    computationTest.augmentMetapathwayNoComputeInverse(types);
    computationTest.addEdges(virtualInputEdges,virtualInputEdgesValues);
    computationTest.addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    computationTest.setDissipationModel(dms2);
    computationTest.setConservationModel(cms2);
    std::vector<double> saturationVector = std::vector<double>(8,0.5);
    std::vector<double> qVector = std::vector<double>(8,0.5);
    //std::vector<double> saturationVector{2,1.5,2,4,1,1,0.5,1.5};
    std::vector<double> result = computationTest.computeAugmentedPerturbationEnhanced2(0,false,std::vector<double>(),qVector);
    std::vector<double> expected{3.256667/2-(input[0]/8),2.736364/2-(input[1]/8),2.924242/2-(input[2]/8),5.593939/2-(input[3]/8),0,0,1.303030/2,2.796970/2};
    // std::vector<double> expected{3.256667/2-(input[0]/4),2.736364/2-(input[1]/4),2.924242/2-(input[2]/4),5.593939/2-(input[3]/4),0,0,1.303030/2,2.796970/2};
    ASSERT_EQ(result.size(),expected.size());
    for (uint i = 0; i < expected.size() ; i++) {
        EXPECT_NEAR(result[i],expected[i],1e-4);
    }
}

