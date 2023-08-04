#pragma  once

#include <array>
#include <cstddef>
#include <iostream>
#include <map>
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
#include<cmath>

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

void printVector(std::vector<int> vec);
void printVector(std::vector<double> vec);

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

/**
scale the hyperbolic tangent function, the return value is always < c , the function is also scaled to grow linearly before reaching the transient
*/
double hyperbolicTangentScaled(double xInput, double scaleFactor );

template<typename T>
std::vector<T> arrayToVector(T* array, int size){
    return std::vector<T>(array, array + size);
}

bool file_exists (const std::string& name);
bool fileExistsPath(const std::string& filePath);
bool folderExists(const std::string& folderPath);
std::vector<std::string> splitString(std::string toSplit , std::string delimiter);

std::pair<std::vector<int>,std::vector<std::tuple<int,int,double>>> edgesFileToEdgesListByIndex(std::string filename);
std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>> edgesFileToEdgesListAndNodesByName(std::string filename);
std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeMatrixToCellVectors(std::string filename, const std::vector<std::string>& finalNames,bool useEntrez=false);
std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeMatrixToCellVectors(std::string filename, const std::vector<std::string>& finalNames,std::vector<std::string> subTypes= std::vector<std::string>(), bool useEntrez=false);
std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> cellInteractionFileToEdgesListAndNodesByName(std::string filename,bool useEntrez=false);
std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> cellInteractionFileToEdgesListAndNodesByName(std::string filename, std::vector<std::string> subtypes,bool useEntrez=false);
std::vector<double> saturationFileToVector(std::string filename,const std::map<std::string, int>& ensembleToIndexMap);
/**
 * \brief   Return the types taken from the file names in a folder with the extension .tsv
 *          that is if the folder contains the files: A.tsv, B.tsv, C.tsv, D.tsv, E.tsv
 *         the function will return the vector {"A","B","C","D","E"}
 */
std::vector<std::string> getTypesFromFolderFileNames(std::string folderPath);
/**
 * \brief   Return the types taken from the first line of a file
 *          that is if the first line contains: name, A, B, C, D, E
 *         the function will return the vector {"A","B","C","D","E"}
 */
std::vector<std::string> getTypesFromMatrixFile(std::string matrixFilepath);

template<typename T>
std::vector<T> getVectorFromFile(std::string filename){
    std::vector<T> retVec;
    std::ifstream inFile(filename);
    std::string line;
    while (std::getline(inFile, line))
    {
        std::stringstream ss(line);
        T temp;
        ss >> temp;
        retVec.push_back(temp);
    }
    inFile.close();
    return retVec;
}

std::map<std::string, std::string>getEnsembletoEntrezidMap();
std::map<std::string, std::vector<std::string>> getFullNodesDescription();
/**
 * \brief   Return the filenames of all files that have the specified extension
 *          in the specified directory and all subdirectories.
 */
std::vector<std::string> get_all(std::string const & root, std::string const & ext);

void saveNodeValues(std::string folderName,int iteration, std::string cellName, std::vector<double> nodeValues,std::vector<std::string> nodeNames, bool useEntrez=false);

template<typename T>
std::vector<T> vectorNormalization(std::vector<T> vec){
    T norm=0;
    for (int i = 0; i < vec.size(); ++i) {
        norm+=vec[i]*vec[i];
    }
    norm=sqrt(norm);
    for (int i = 0; i < vec.size(); ++i) {
        vec[i]=vec[i]/norm;
    }
    return vec;
}

double vectorNorm(std::vector<double> vec);

template<typename T>
std::vector<T> vectorScalarMultiplication(std::vector<T> vec, T scalar){
    for (uint i = 0; i < vec.size(); ++i) {
        vec[i]=vec[i]*scalar;
    }
    return vec;
}

template<typename T>
bool vectorContains(std::vector<T> vec, T element){
    for (uint i = 0; i < vec.size(); ++i) {
        if(vec[i]==element){
            return true;
        }
    }
    return false;
}