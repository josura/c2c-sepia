#pragma  once

#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <bitset>
#include <random>
#include <chrono>
#include <algorithm>

#include "WeightedEdgeGraph.h"

#define INTMAX 100
#define DOUBLEMAX 100.0

//custom types
typedef std::vector<int> NodeList;
using NodeBitList = std::vector<bool>;
using NodeBitArray = bool*;
using NodeSet = std::unordered_set<int>;


std::ostream& operator<< (std::ostream &out, WeightedEdgeGraph const& data);  

std::ostream& operator<< (std::ostream &out, NodeBitList const& data);
std::ostream& operator<< (std::ostream &out, NodeList const& data);
std::ostream& operator<< (std::ostream &out, NodeSet const& data);

NodeList* nodeBitArrayToList(NodeBitArray const& nodeArray,int arraySize);
NodeSet* nodeBitArrayToSet(NodeBitArray const& nodeArray,int arraySize);

void printNodeBitArray(NodeBitArray nodeArray,int size);

int randomNumber(int min, int max);
double randomRealNumber(double min, double max);
char generateRandomCharacter();

std::vector<int> randomVector(int min, int max , int size);

NodeBitArray randomBooleanArray(int size);

void printUsage(std::string execName);

std::string nodeBitArrayToString(NodeBitArray nodeArray,int size);

std::unordered_set<int> intersectionSet(std::unordered_set<int> set1,std::unordered_set<int> set2);

void setRandom(int& val);
void setRandom(double& val);
void setRandom(char& val);