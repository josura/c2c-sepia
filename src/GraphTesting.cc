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

TEST_F(GraphTesting, assignmentWorks2) {
  WeightedEdgeGraph graphTest;
  graphTest.assign(*g4_);  //assignment is done in this way since the = operator is working only on dereferenced objects, that is *obj1=*obj2

  ASSERT_EQ(graphTest.getNumNodes(), 5);
  ASSERT_EQ(graphTest.getNumEdges(), 2);

  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight("node1", "node3"), 2.3);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight("node3", "node2"), 0.3);


  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(graphTest.getIndexFromName("node1") , graphTest.getIndexFromName("node3")), 2.3);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(graphTest.getIndexFromName("node3"), graphTest.getIndexFromName("node2")), 0.3);

}

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
  WeightedEdgeGraph graphTest;
  graphTest.assign(*g5_);

  graphTest.addEdge(1,2,0.1)->addEdge(1,3,0.2)->addEdge(3,0,0.4);

  EXPECT_EQ(graphTest.getNumEdges(), 3);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(1, 2), 0.1);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(1, 3), 0.2);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(3, 0), 0.4);
}

TEST_F(GraphTesting, addingEdgesNames){
  WeightedEdgeGraph graphTest;
  graphTest.assign(*g5_);

  graphTest.addEdge("node2","node3",0.1)->addEdge("node2","node4",0.2)->addEdge("node4","node1",0.4);

  EXPECT_EQ(graphTest.getNumEdges(), 3);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight("node2", "node3"), 0.1);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight("node2", "node4"), 0.2);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight("node4", "node1"), 0.4);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(graphTest.getIndexFromName("node2"),graphTest.getIndexFromName("node3") ), 0.1);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(graphTest.getIndexFromName("node2"), graphTest.getIndexFromName("node4")), 0.2);
  EXPECT_FLOAT_EQ(graphTest.getEdgeWeight(graphTest.getIndexFromName("node4"), graphTest.getIndexFromName("node1")), 0.4);

}

TEST_F(GraphTesting, addingNodeOnlyValueAndDefualt){
  WeightedEdgeGraph graphTest;
  graphTest.assign(*g5_);

  graphTest.addEdge("node2","node3",0.1)->addEdge("node2","node4",0.2)->addEdge("node4","node1",0.4);
  graphTest.addNode(4.0)->addNode();
  int nodeDefaultIndex = graphTest.getNumNodes()-1;
  int nodeValueIndex = nodeDefaultIndex-1;
  std::string nodeValueName = std::to_string(nodeValueIndex);
  std::string nodeDefaultName =  std::to_string(nodeDefaultIndex);
  ASSERT_EQ(graphTest.getNumEdges(), 3);
  ASSERT_EQ(graphTest.getNumNodes(), 7);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeValueIndex),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeDefaultIndex),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeValueName),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeDefaultName),0.0);
  EXPECT_NO_THROW({graphTest.getAdjList(nodeValueIndex);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeDefaultIndex);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeValueName);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeDefaultName);});
  EXPECT_EQ(graphTest.adjMatrix.getCols(), graphTest.getNumNodes());
  EXPECT_EQ(graphTest.adjMatrix.getRows(), graphTest.getNumNodes());

}

TEST_F(GraphTesting, addingNodeNameAndValue){
  WeightedEdgeGraph graphTest;
  graphTest.assign(*g5_);
  std::string nodeValueName = "nodetest";
  std::string nodeDefaultName =  "nodetest2";

  graphTest.addEdge("node2","node3",0.1)->addEdge("node2","node4",0.2)->addEdge("node4","node1",0.4);
  graphTest.addNode("nodetest",4.0)->addNode("nodetest2");
  
  int nodeValueIndex = graphTest.getIndexFromName(nodeValueName);
  int nodeDefaultIndex = graphTest.getIndexFromName(nodeDefaultName);
  ASSERT_EQ(graphTest.getNumEdges(), 3);
  ASSERT_EQ(graphTest.getNumNodes(), 7);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeValueIndex),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeDefaultIndex),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeValueName),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeDefaultName),0.0);
  EXPECT_NO_THROW({graphTest.getAdjList(nodeValueIndex);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeDefaultIndex);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeValueName);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeDefaultName);});
  EXPECT_EQ(graphTest.adjMatrix.getCols(), graphTest.getNumNodes());
  EXPECT_EQ(graphTest.adjMatrix.getRows(), graphTest.getNumNodes());
}


TEST_F(GraphTesting, addingNodesOnlyValues){
  WeightedEdgeGraph graphTest;
  graphTest.assign(*g5_);
  std::vector<double> values{4.0,8.0,3.2};
  graphTest.addEdge("node2","node3",0.1)->addEdge("node2","node4",0.2)->addEdge("node4","node1",0.4);
  graphTest.addNodes(values);
  int nodesStartIndex = graphTest.getNumNodes()-3;
  std::vector<std::string> nodeStartName{std::to_string(nodesStartIndex),std::to_string(nodesStartIndex+1),std::to_string(nodesStartIndex+2)};
  ASSERT_EQ(graphTest.getNumEdges(), 3);
  ASSERT_EQ(graphTest.getNumNodes(), 8);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+1),8.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+2),3.2);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeStartName[0]),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeStartName[1]),8.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodeStartName[2]),3.2);
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex+1);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex+2);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeStartName[0]);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeStartName[1]);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodeStartName[2]);});
  EXPECT_EQ(graphTest.adjMatrix.getCols(), graphTest.getNumNodes());
  EXPECT_EQ(graphTest.adjMatrix.getRows(), graphTest.getNumNodes());

  std::vector<double> graphNodeValues = graphTest.getNodeValues(nodeStartName);
  for (uint i = 0; i < nodeStartName.size(); i++) {
    EXPECT_FLOAT_EQ(graphNodeValues[i], values[i]) << "Vectors nodevalues and graphtest differ at index " << i << " at node:" << nodeStartName[i];
  }

}

TEST_F(GraphTesting, addingNodesNamesAndValuesAndDefault){
  WeightedEdgeGraph graphTest;
  graphTest.assign(*g5_);
  std::vector<double> values{4.0,8.0,3.2};
  std::vector<std::string> namesWithVal{"nodetest1","nodetest2","nodetest3"};
  std::vector<std::string> namesWithDefault{"nodetest4","nodetest5","nodetest6","nodetest7"};

  graphTest.addEdge("node2","node3",0.1)->addEdge("node2","node4",0.2)->addEdge("node4","node1",0.4);
  graphTest.addNodes(namesWithVal,values);


  int nodesStartIndex = graphTest.getNumNodes()-3;

  ASSERT_EQ(graphTest.getNumNodes(), 8);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+1),8.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+2),3.2);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithVal[0]),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithVal[1]),8.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithVal[2]),3.2);
  EXPECT_EQ(graphTest.adjMatrix.getCols(), graphTest.getNumNodes());
  EXPECT_EQ(graphTest.adjMatrix.getRows(), graphTest.getNumNodes());

  graphTest.addNodes(namesWithDefault);
  nodesStartIndex = graphTest.getNumNodes()-7;

  
  ASSERT_EQ(graphTest.getNumNodes(), 12);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+1),8.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+2),3.2);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+3),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+4),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+5),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(nodesStartIndex+6),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithVal[0]),4.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithVal[1]),8.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithVal[2]),3.2);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithDefault[0]),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithDefault[1]),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithDefault[2]),0.0);
  EXPECT_FLOAT_EQ(graphTest.getNodeValue(namesWithDefault[3]),0.0);
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex+1);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex+2);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex+3);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex+4);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex+5);});
  EXPECT_NO_THROW({graphTest.getAdjList(nodesStartIndex+6);});
  EXPECT_NO_THROW({graphTest.getAdjList(namesWithVal[0]);});
  EXPECT_NO_THROW({graphTest.getAdjList(namesWithVal[1]);});
  EXPECT_NO_THROW({graphTest.getAdjList(namesWithVal[2]);});
  EXPECT_NO_THROW({graphTest.getAdjList(namesWithDefault[0]);});
  EXPECT_NO_THROW({graphTest.getAdjList(namesWithDefault[1]);});
  EXPECT_NO_THROW({graphTest.getAdjList(namesWithDefault[2]);});
  EXPECT_NO_THROW({graphTest.getAdjList(namesWithDefault[3]);});
  EXPECT_EQ(graphTest.adjMatrix.getCols(), graphTest.getNumNodes());
  EXPECT_EQ(graphTest.adjMatrix.getRows(), graphTest.getNumNodes());
  //also test if getNodeValues with default argument returns all the nodes with the right values
  std::vector<double> graphNodeValues = graphTest.getNodeValues();
  for (uint i = 0; i < namesWithVal.size(); i++) {
    EXPECT_FLOAT_EQ(graphNodeValues[nodesStartIndex + i], values[i]) << "Vectors nodevalues and graphtest differ at index " << i << " at node:" << namesWithVal[i];
  }
  for (int i = 0; i < nodesStartIndex; i++) {
    EXPECT_FLOAT_EQ(graphNodeValues[i],nodeValues[i]);
  }
}

TEST_F(GraphTesting, setNodeValueIndex){
  g5_->setNodeValue(1,0.2);
  EXPECT_FLOAT_EQ(g5_->getNodeValue(1),0.2);
}
TEST_F(GraphTesting, setNodeValueName){
  g5_->setNodeValue("node2",4.5);
  EXPECT_FLOAT_EQ(g5_->getNodeValue("node2"),4.5);
}
//throws and unexpected behaviour management TODO

TEST_F(GraphTesting, gettingNodeValueOfNotPresentNode){
  EXPECT_ANY_THROW(g5_->getNodeValue("nodenotPresent"));
  EXPECT_ANY_THROW(g5_->getNodeValue(g5_->getNumNodes()));
  EXPECT_ANY_THROW(g5_->getNodeValue(-1));
}
TEST_F(GraphTesting, gettingNodesValuesOfAtLeastOneNotPresentNode){
  EXPECT_ANY_THROW(g5_->getNodeValues(std::vector<std::string>{"nodeNotInpis","node1","node2"}));
  EXPECT_NO_THROW(g5_->getNodeValues(std::vector<std::string>{"node1","node2"}));
}
TEST_F(GraphTesting, addingEdgeOfNotPresentNode){
  EXPECT_ANY_THROW(g5_->addEdge("node1","nodeNot",3.2));
}
TEST_F(GraphTesting, addingEdgeOfNotPresentNodes){
  EXPECT_ANY_THROW(g5_->addEdge("nodenot1","nodeNot",4.2));
}
TEST_F(GraphTesting, gettingAdjListOfNotPresentNode){
  EXPECT_ANY_THROW(g5_->getAdjList(6));
  EXPECT_ANY_THROW(g5_->getAdjList(5));
  EXPECT_ANY_THROW(g5_->getAdjList(-1));
  EXPECT_ANY_THROW(g5_->getAdjList("nodenot"));
}

TEST_F(GraphTesting, setNodeValueOfNotPresentNode){
  EXPECT_ANY_THROW(g5_->setNodeValue(-1,0.2));
  EXPECT_ANY_THROW(g5_->setNodeValue(6,0.2));
  EXPECT_ANY_THROW(g5_->setNodeValue(5,0.2));
}