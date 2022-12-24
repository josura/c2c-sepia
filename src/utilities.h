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

//custom types
typedef std::vector<uint> NodeList;
using NodeBitList = std::vector<bool>;
using NodeBitArray = bool*;
using NodeSet = std::unordered_set<uint>;


std::ostream& operator<< (std::ostream &out, WeightedEdgeGraph const& data);

std::ostream& operator<< (std::ostream &out, NodeBitList const& data);
std::ostream& operator<< (std::ostream &out, NodeList const& data);
std::ostream& operator<< (std::ostream &out, NodeSet const& data);

NodeList* nodeBitArrayToList(NodeBitArray const& nodeArray,uint arraySize);
NodeSet* nodeBitArrayToSet(NodeBitArray const& nodeArray,uint arraySize);

void printNodeBitArray(NodeBitArray nodeArray,uint size);

int randomNumber(int min, int max);
double randomRealNumber(double min, double max);

std::vector<uint> randomVector(int min, int max , uint size);

NodeBitArray randomBooleanArray(uint size);

void printUsage(std::string execName);

std::string nodeBitArrayToString(NodeBitArray nodeArray,uint size);

std::unordered_set<uint> intersectionSet(std::unordered_set<uint> set1,std::unordered_set<uint> set2);