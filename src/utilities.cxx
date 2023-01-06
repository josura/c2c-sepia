#include "utilities.h"
#include "WeightedEdgeGraph.h"
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <random>
#include <tuple>
#include <unordered_set>
#include <vector>

using namespace std;

std::ostream& operator<< (std::ostream &out, WeightedEdgeGraph const& data) {
            out << data.getNumNodes() << " " << data.getNumEdges() <<std::endl;
            string nodeweights = data.getNodeWeightsStr();
            out << nodeweights << std::endl;
            out << "Adj Lists" << std::endl;
            for(int i = 0; i<data.getNumNodes(); i++){
                out << "node " << i << " :" << data.getAdjListStr(i) << std::endl;
            }
            out << "Edges vector: {";
            for(auto it = data.getEdgesVector()->cbegin();it!=data.getEdgesVector()->cend();it++){
                out << "(" << get<0>(*it)<< ","<< get<1>(*it)<< ","<< get<2>(*it)<< "," << ")," ;
            }
            out << "}"<< std::endl;
            return out;
        }



//random generation for different types

int randomNumber(int min, int max){
    // fixed seed because repeatability
    //unsigned seed = 777;
    std::random_device r;
    // range [min,max[
    std::default_random_engine e1(r());
    max = (max-1<0) ? 0 : max-1;
    std::uniform_int_distribution<int> uniform_dist(min, max);
    return uniform_dist(e1);
}

double randomRealNumber(double min, double max){
    // fixed seed because repeatability
    //unsigned seed = 777;
    std::random_device r;
    // range [min,max]
    std::default_random_engine e1(r());
    std::uniform_real_distribution<double> uniform_dist(min, max);
    return uniform_dist(e1);
}

char generateRandomCharacter() {
  // Generate a random number between 0 and 25
  int randomNumberInt = randomNumber(0, 25);

  // Convert the random number to a character
  // by adding 'a' (the ASCII value for 'a' is 97)
  char randomCharacter = (char)(randomNumberInt + 'a');

  return randomCharacter;
}

void createRandom(int& val) { 
    val =  randomNumber(-INTMAX, INTMAX);
}
void createRandom(double& val) { 
    val = randomRealNumber(-DOUBLEMAX, DOUBLEMAX);
}
void createRandom(char& val) { 
    val = generateRandomCharacter();
}