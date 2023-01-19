#include <gtest/gtest.h>
#include <string>
#include <vector>
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

    g4_ = new WeightedEdgeGraph(nodeNames);
    g5_ = new WeightedEdgeGraph(nodeNames,nodeValues);

    g4_->addEdge("node1","node3",2.3)->addEdge("node3","node2",0.3);
     
  }
  void TearDown() override{
    delete g0_;
    delete g1_;
    delete g2_;
    delete g3_;
    delete g4_;
    delete g5_;
  }

  // void TearDown() override {}

  WeightedEdgeGraph* g0_;       //testing default constructor
  WeightedEdgeGraph* g1_;    //testing single parameter constructor
  WeightedEdgeGraph* g2_;  //testing passing another weighted edge graph as input to constructor (= operator)
  WeightedEdgeGraph* g3_;  //testing passing an adjacency matrix
  WeightedEdgeGraph* g4_;  //testing passing a vector for the node names
  WeightedEdgeGraph* g5_;  //testing passing a vector for the node names and a vector for the node values


  std::vector<std::string> nodeNames{"node1","node2","node3","node4","node5"};
  std::vector<double> nodeValues{0.3,4.1,3.8,8.2,9.5};
  
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

TEST_F(GraphTesting, assignmentWorks1) {
  EXPECT_EQ(g2_->getNumNodes(), 4);
  EXPECT_EQ(g2_->getNumEdges(), 2);
}

// TEST_F(GraphTesting, assignmentWorks2) {
//   WeightedEdgeGraph graphTest = *g4_;

//   ASSERT_EQ(graphTest.getNumNodes(), 5);
//   ASSERT_EQ(graphTest.getNumEdges(), 2);

//   EXPECT_FLOAT_EQ(graphTest.getEdgeWeight("node1", "node3"), 2.3);
//   EXPECT_FLOAT_EQ(graphTest.getEdgeWeight("node3", "node2"), 0.3);


//   EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(graphTest.getIndexFromName("node1") , graphTest.getIndexFromName("node3")), 2.3);
//   EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(graphTest.getIndexFromName("node3"), graphTest.getIndexFromName("node2")), 0.3);

// }

TEST_F(GraphTesting, assignmentAdjMatrixWorks) {
  EXPECT_EQ(g3_->getNumNodes(), 4);
  EXPECT_EQ(g3_->getNumEdges(), 2);
}


TEST_F(GraphTesting, constructorNodeNames) {
  EXPECT_EQ(g4_->getNumNodes(), 5);
  EXPECT_EQ(g4_->getNumEdges(), 2);
}


TEST_F(GraphTesting, constructorNodeNamesAndValues) {
  EXPECT_EQ(g5_->getNumNodes(), 5);
  EXPECT_EQ(g5_->getNumEdges(), 0);
  ASSERT_EQ(nodeNames.size(), nodeValues.size()) << "node names and node values have unequal length";
  auto g5Map = g5_->getNodeToIndexMap();
  ASSERT_EQ(nodeNames.size(), g5Map.size()) << "node names and nodes in graph have unequal length";

  for (uint i = 0; i < nodeNames.size(); ++i) {
    EXPECT_FLOAT_EQ(nodeValues[i], g5_->getNodeValue(nodeNames[i])) << "Vectors nodevalues and g5Map differ at index " << i;
  }

}

TEST_F(GraphTesting, addingEdgesIndex){
  // WeightedEdgeGraph graphTest = *g5_;

  // graphTest.addEdge(1,2,0.1)->addEdge(1,3,0.2)->addEdge(3,0,0.4);

  // EXPECT_EQ(graphTest.getNumEdges(), 3);
  // EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(1, 2), 0.1);
  // EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(1, 3), 0.2);
  // EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(3, 0), 0.4);
}

// TEST_F(GraphTesting, addingEdgesNames){

// }
