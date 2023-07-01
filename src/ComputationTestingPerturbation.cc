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
            c0  = new Computation();
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
        std::vector<double> virtualInputEdgesValues = {0.4,0.5,0.7,0.2};
        std::vector<std::pair<std::string, std::string>> virtualOutputEdges{{"node2","v-out:type2"},
                                                                            {"node4","v-out:type3"}
                                                                                     };
        std::vector<double> virtualOutputEdgesValues = {0.4,0.4};

        Computation* c0;       //testing default constructor
        Computation* c1;       //testing general constructor

        DissipationModel* dms = new DissipationModelScaled( [](double x){return 0;} );
        DissipationModel* dms = new DissipationModelScaled(); //default constructor uses the functions that returns 0.5 for every iteration

        ConservationModel* cms = new ConservationModel( [](double x){return 0;} );
        ConservationModel* cms = new ConservationModel(); //default constructor uses the functions that returns 0.5 for every iteration
  
};

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationNoConservationDefault) {

}


TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationNoConservation) {

}


TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectDissipationDefaultNoConservation) {

}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectDissipationNoConservation) {

}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationConservationDefaultQempty) {

}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationConservationQempty) {

}

TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectNoDissipationConservationQpassed) {

}


TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectDissipationConservationQempty) {

}


TEST_F(ComputationTestingPerturbation, computePerturbationIsCorrectDissipationConservationQpassed) {

}

