#include <gtest/gtest.h>
#include "PropagationModel.hxx"
#include "PropagationModelNeighbors.hxx"
#include "WeightedEdgeGraph.h"
#include <armadillo>
#include <functional>

class PropagationModelTesting : public ::testing::Test {
 protected:
  void SetUp() override {

    input_[0] = 1;
    input_[4] = 1;

    q0_ = new WeightedEdgeGraph(6);
    q1_ = new WeightedEdgeGraph(6);
    q2_ = new WeightedEdgeGraph(6);


    q1_->addEdge(0,1,1);
    q1_->addEdge(1,2,1);
    q1_->addEdge(2,3,1);
    q1_->addEdge(3,4,1);
    q1_->addEdge(4,5,1);

    q2_->addEdge(0,1,1);
    q2_->addEdge(1,2,1);
    q2_->addEdge(2,3,1);
    q2_->addEdge(3,4,1);
    q2_->addEdge(4,5,1);

    // q0_ is a graph with 6 nodes, 0 edges
    // q1_ is a graph with 6 nodes, 5 edges
    m0_  = new PropagationModelNeighbors(q0_);
    m1_  = new PropagationModelNeighbors(q1_,[](double time)-> double{return 1;});
    m2_  = new PropagationModelNeighbors(q2_);

    


     
  }
  void TearDown() override{
    delete m0_;
  }

  // void TearDown() override {}
  PropagationModel* m0_;       //testing default constructor for neighbors model
  PropagationModel* m1_;       //testing constructor with scale function for neighbors model
  PropagationModel* m2_;       //testing default constructor with edges for neighbors model

  WeightedEdgeGraph *q0_;   // using graph with 6 nodes
  WeightedEdgeGraph *q1_;   
  WeightedEdgeGraph *q2_;   // using graph with 6 nodes and 5 edges

  arma::Col<double> input_{6,arma::fill::zeros};

  
};

TEST_F(PropagationModelTesting, propagateWorksWithDefaultScaleFunction) {
  arma::Col<double> output = m0_->propagate(input_,0);
  EXPECT_DOUBLE_EQ(output(0), 1);
  EXPECT_DOUBLE_EQ(output(1), 0);
  EXPECT_DOUBLE_EQ(output(2), 0);
  EXPECT_DOUBLE_EQ(output(3), 0);
  EXPECT_DOUBLE_EQ(output(4), 1);
  EXPECT_DOUBLE_EQ(output(5), 0);
}

TEST_F(PropagationModelTesting, propagateWorksWithDefaultScaleFunctionWithEdges) {
  arma::Col<double> output = m2_->propagate(input_,0);
  std::function<double(double)> scaleFunc = [](double time)-> double{return 0.5;};
  EXPECT_DOUBLE_EQ(output(0), 1);
  EXPECT_DOUBLE_EQ(output(1), scaleFunc(1));
  EXPECT_DOUBLE_EQ(output(2), 0);
  EXPECT_DOUBLE_EQ(output(3), 0);
  EXPECT_DOUBLE_EQ(output(4), 1);
  EXPECT_DOUBLE_EQ(output(5), scaleFunc(1));
}

TEST_F(PropagationModelTesting, propagateWorksWithScaleFunction) {
  arma::Col<double> output = m1_->propagate(input_,0);
  EXPECT_DOUBLE_EQ(output(0), 1);
  EXPECT_DOUBLE_EQ(output(1), 1);
  EXPECT_DOUBLE_EQ(output(2), 0);
  EXPECT_DOUBLE_EQ(output(3), 0);
  EXPECT_DOUBLE_EQ(output(4), 1);
  EXPECT_DOUBLE_EQ(output(5), 1);
}