#include <gtest/gtest.h>
#include "WeightedEdgeGraph.h"

class GraphTesting : public ::testing::Test {
 protected:
  void SetUp() override {
    // q0_ remains empty
    g0_  = new WeightedEdgeGraph();
    g1_ = new WeightedEdgeGraph(4);
    g2_ = new WeightedEdgeGraph();

    
    g1_->addEdge(1,2,0.3)->addEdge(0,2,0.4);
    *g2_ = *g1_;   //Problem with the equal operator(assignment)

    g3_ = new WeightedEdgeGraph(g1_->adjMatrix);
     
  }
  void TearDown() override{
    delete g0_;
    delete g1_;
    delete g2_;
    delete g3_;
  }

  // void TearDown() override {}

  WeightedEdgeGraph* g0_;       //testing default constructor
  WeightedEdgeGraph* g1_;    //testing single parameter constructor
  WeightedEdgeGraph* g2_;  //testing passing another weighted edge graph as input to constructor (= operator)
  WeightedEdgeGraph* g3_;  //testing passing an adjacency matrix
};

TEST_F(GraphTesting, constructorWorksDefault) {
  EXPECT_EQ(g0_->getNumNodes(), 0);
  EXPECT_EQ(g0_->getNumEdges(), 0);
}


TEST_F(GraphTesting, constructorWorksNumNodes) {
  EXPECT_EQ(g1_->getNumNodes(), 4);
  EXPECT_EQ(g1_->getNumEdges(), 2);
}


TEST_F(GraphTesting, constructorWorksNumNodesSquareMatrix) {
  EXPECT_EQ(g1_->adjMatrix.getCols(), g1_->adjMatrix.getRows());
}


TEST_F(GraphTesting, constructorWorksMatrixHasRightEdges) {
  EXPECT_FLOAT_EQ(g1_->adjMatrix.getValue(1, 2),0.3);
  EXPECT_FLOAT_EQ(g1_->adjMatrix.getValue(0, 2),0.4);
}

TEST_F(GraphTesting, assignmentWorks) {
  EXPECT_EQ(g2_->getNumNodes(), 4);
  EXPECT_EQ(g2_->getNumEdges(), 2);
}

TEST_F(GraphTesting, assignmentAdjMatrixWorks) {
  EXPECT_EQ(g3_->getNumNodes(), 4);
  EXPECT_EQ(g3_->getNumEdges(), 2);
}