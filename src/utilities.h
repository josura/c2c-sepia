#pragma  once

#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>
#include <bitset>
#include <random>
#include <chrono>
#include <algorithm>
#include <limits>    // for std::numeric_limits
#include <stdexcept> // for std::overflow_error
#include <sys/stat.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


#define INTMAX 100
#define DOUBLEMAX 100.0

//custom types
typedef std::vector<int> NodeList;
using NodeBitList = std::vector<bool>;
using NodeBitArray = bool*;
using NodeSet = std::unordered_set<int>; 

std::ostream& operator<< (std::ostream &out, NodeBitList const& data);
std::ostream& operator<< (std::ostream &out, NodeList const& data);
std::ostream& operator<< (std::ostream &out, NodeSet const& data);

NodeList* nodeBitArrayToList(NodeBitArray const& nodeArray,int arraySize);
NodeSet* nodeBitArrayToSet(NodeBitArray const& nodeArray,int arraySize);

void printNodeBitArray(NodeBitArray nodeArray,int size);

//random generation for integers
int randomNumber(int min, int max);
//random generation for doubles
double randomRealNumber(double min, double max);
//random generation for character
char generateRandomCharacter();

//generate random vector
std::vector<int> randomVector(int min, int max , int size);

NodeBitArray randomBooleanArray(int size);

void printUsage(std::string execName);

std::string nodeBitArrayToString(NodeBitArray nodeArray,int size);

std::unordered_set<int> intersectionSet(std::unordered_set<int> set1,std::unordered_set<int> set2);

void setRandom(int& val);
void setRandom(double& val);
void setRandom(char& val);

//generate random matrix
/*
@param rows: number of rows
@param cols: number of cols
*/
template<typename T>
std::vector<std::vector< T>> createRandomMatrixVector(int rows,int cols){
    std::vector<std::vector< T>> retMat=std::vector<std::vector<double>>(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            setRandom(retMat[i][j]);
        }
    }
    return retMat;
}

template<typename T>
int getIndex(std::vector<T> v, T K)
{
    auto it = find(v.begin(), v.end(), K);
  
    // If element was found
    if (it != v.end()) 
    {
        return (it - v.begin());
    }
    else {
        return -1;
    }
}

// convert size to int and  launch an exception if it is not possible
int SizeToInt(size_t u);

bool approximatelyEqual(float a, float b, float epsilon);
bool essentiallyEqual(float a, float b, float epsilon);
bool definitelyGreaterThan(float a, float b, float epsilon);
bool definitelyLessThan(float a, float b, float epsilon);


bool approximatelyEqual(double a, double b, double epsilon);
bool essentiallyEqual(double a, double b, double epsilon);
bool definitelyGreaterThan(double a, double b, double epsilon);
bool definitelyLessThan(double a, double b, double epsilon);

template<typename T>
std::vector<T> arrayToVector(T* array, int size){
    return std::vector<T>(array, array + size);
}

inline bool file_exists (const std::string& name);
std::vector<std::string> splitString(std::string toSplit , std::string delimiter);

std::pair<std::vector<int>,std::vector<std::tuple<int,int,double>>> edgesFileToEdgesListByIndex(std::string filename);
std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>> edgesFileToEdgesListAndNodesByName(std::string filename);