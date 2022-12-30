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
    g2_ = g1_;
     
  }
  void TearDown() override{
    delete g0_;
    delete g1_;
    delete g2_;
  }

  // void TearDown() override {}

  WeightedEdgeGraph* g0_;       //testing default constructor
  WeightedEdgeGraph* g1_;    //testing single parameter constructor
  WeightedEdgeGraph* g2_;  //testing passing another weighted edge graph as input to constructor (= operator)
};

TEST_F(GraphTesting, constructorWorksDefault) {
  EXPECT_EQ(g0_->getNumNodes(), 0);
  EXPECT_EQ(g0_->getNumEdges(), 0);
}


TEST_F(GraphTesting, constructorWorksNumNodes) {
  EXPECT_EQ(g1_->getNumNodes(), 4);
  EXPECT_EQ(g1_->getNumEdges(), 2);
}
