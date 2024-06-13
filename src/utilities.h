#pragma  once

#include <array>
#include <cstddef>
#include <exception>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
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

// hash function for unordered_map
#include <bits/stdc++.h>

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

template<typename T>
void printVector(std::vector<T> vec){
    for (uint i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << " , ";
    }
    std::cout << std::endl;
}

/**
 * \brief   Generate a random integer number between min and max
 * \return  the random number
*/
int randomNumber(int min, int max);
/**
 * \brief   Generate a random real number between min and max
 * \return  the random number
*/
double randomRealNumber(double min, double max);
/**
 * \brief   Generate a random character between A and Z
 * \return  the random character
*/
char generateRandomCharacter();

/**
 * \brief   Generate a random integer vector of length size
 * \return  the random vector
*/
std::vector<int> randomVector(int min, int max , int size);

/**
 * \brief   Generate a random Boolean vector of length size
 * \return  the random vector
*/
NodeBitArray randomBooleanArray(int size);


long int szudzik(int x, int y);

void printUsage(std::string execName);

std::string nodeBitArrayToString(NodeBitArray nodeArray,int size);

/**
 * \brief  Generate the intersection of two sets
 * \return  the intersection of the two sets
*/
std::unordered_set<int> intersectionSet(std::unordered_set<int> set1,std::unordered_set<int> set2);

void setRandom(int& val);
void setRandom(double& val);
void setRandom(char& val);

//generate random matrix
/**
 * \brief   Generate a random matrix of size rows x cols
 * \return  the random matrix
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

template<typename T>
bool controlForDuplicates(std::vector<T> v){
    std::vector<T> v2 = v;
    std::sort(v2.begin(), v2.end());
    auto last = std::unique(v2.begin(), v2.end());
    v2.erase(last, v2.end());
    if(v.size()!=v2.size()){
        return true;
    }
    return false;
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

/**
 * \brief  Linear interpolation between two values, a and b, with a parameter t going from 0 to 1
 * \return  the interpolated value
*/

double lerping(double a, double b, double t);

/**
 * \brief  Convert an array to a vector
 * \return  the vector
*/

template<typename T>
std::vector<T> arrayToVector(T* array, int size){
    return std::vector<T>(array, array + size);
}

/**
 * \brief   Control if the file exists in the current directory
 * \return  true if the file exists, false otherwise
 */
bool file_exists (const std::string& name);

/**
 * \brief   Control if the file exists in the specified path
 * \return  true if the file exists, false otherwise
 */
bool fileExistsPath(const std::string& filePath);
/**
 * \brief   Control if the folder exists in the specified path
 * \return  true if the folder exists, false otherwise
 */
bool folderExists(const std::string& folderPath);
/**
 * \brief   Create the folder in the specified path
 * \return  true if the folder is created, false otherwise
 */
bool createFolder(const std::string& folderPath);
/**
 * \brief  Split a string into a vector of strings using the specified delimiter
 * \return  the vector of strings
*/
std::vector<std::string> splitStringIntoVector(std::string toSplit , std::string delimiter);

/**
 * \brief  Split a string into a vector of strings using the specified delimiter into two part, the first part is the string before the delimiter and the second part is the string after the delimiter
 * the first delimiter found is used
 * \return  the vector of strings
*/
std::vector<std::string> splitStringIntoVectorTwoParts(std::string toSplit , std::string delimiter);

/**
 * \brief  Split a string representing a virtual node (v-in or v-out) into a vector of strings, where the first element is the type and the second element is the name of the node if present
 * \return  the vector of strings
*/
std::vector<std::string> splitVirtualNodeStringIntoVector(std::string toSplit);

std::pair<std::vector<int>,std::vector<std::tuple<int,int,double>>> edgesFileToEdgesListByIndex(std::string filename);
std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>> edgesFileToEdgesListAndNodesByName(std::string filename);
std::pair<std::vector<std::string>,std::vector<std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>>>> edgesFileToEdgesListAndNodesByNameFromFolder(std::string filename);
std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeMatrixToCellVectors(std::string filename, const std::vector<std::string>& finalNames,bool useEntrez=false);
std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeMatrixToCellVectors(std::string filename, const std::vector<std::string>& finalNames,std::vector<std::string> subType, bool useEntrez=false);
std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeCellVectorsFromFolder(std::string folderPath,const std::vector<std::string>& allTypes, const std::vector<std::vector<std::string>>& finalNames,std::vector<std::string> subType = std::vector<std::string>(), bool useEntrez=false);
std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> interactionFileToEdgesListAndNodesByName(std::string filename,bool useEntrez=false);
std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> interactionFileToEdgesListAndNodesByName(std::string filename, std::vector<std::string> subtypes,bool useEntrez=false);

std::map<std::string,std::vector<std::string>> nodeNamesFromFolder(std::string folderPath);

// function to return the keys in a map
template<typename K, typename V>
std::vector<K> getKeys(std::map<K,V> const& input_map) {
    std::vector<K> retval;
    for (auto const& element : input_map) {
        retval.push_back(element.first);
    }
    return retval;
}

// A hash function used to hash a pair of any kind
struct hash_pair_strings {
    size_t operator()(const std::pair<std::string, std::string>& p) const
    {
        std::string tmp;
        // in case the order of the pair does not matter
        // if( p.first < p.second ) {
        //     tmp = p.first + p.second;
        // }
        // else {
        //     tmp = p.second + p.first;
        // }

        // in case the order of the pair matters
        tmp = p.first + p.second;

        auto hashStrings = std::hash<std::string>();
         
        // If hash1 == hash2, their XOR is zero.
        return hashStrings(tmp);
    }
};

struct hash_pair_ints {
    size_t operator()(const std::pair<int, int>& p) const
    {
        long int tmp;
        tmp = szudzik(p.first,p.second);

        auto hashStrings = std::hash<long int>();
         
        // If hash1 == hash2, their XOR is zero.
        return hashStrings(tmp);
    }
};

struct hash_quadruple_strings {
    size_t operator()(const std::tuple<std::string, std::string,std::string, std::string>& t) const
    {
        std::string tmp;
        // in case the order of the pair does not matter
        // if( p.first < p.second ) {
        //     tmp = p.first + p.second;
        // }
        // else {
        //     tmp = p.second + p.first;
        // }

        // in case the order of the tuple matters
        tmp = std::get<0>(t) + std::get<1>(t) + std::get<2>(t) + std::get<3>(t);

        auto hashStrings = std::hash<std::string>();
         
        // If hash1 == hash2, their XOR is zero.
        return hashStrings(tmp);
    }
};

/**
 * \brief   Returns a boolean value if the set contains the interval(width is nonzero):
 *         - set is the set of values
 *         - lower is the lower bound of the interval
 *         - upper  is the upper bound of the interval
 * \return true if the set contains the interval, false otherwise
 */
bool setDoubleContainsInterval(std::set<double> set, double lower, double upper);


/**
 * \brief   Returns how many values fall in the interval:
 *         - set is the set of values
 *         - lower is the lower bound of the interval
 *         - upper  is the upper bound of the interval
 * \return the number of values that fall in the interval
 */
int setDoubleIntervalWidth(std::set<double> set, double lower, double upper);

/**
 * \brief   Returns the new virtual nodes associated with a type along the edges in the augmented graph, it also return the graph of the interactions between types(as a tuple of:
 *         - the start type/agent
 *         - the end type/agent
 *         - the contact times of the interaction
 *         granularity needs to be specified as an argument
 * \return the pair (map of the new virtual nodes associated with a type, graph of the interactions between types)
 */
std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::unordered_set<int>, double>>> interactionContactsFileToEdgesListAndNodesByName(std::string filename, std::vector<std::string> subtypes,int maximumIntertypeTime=INT32_MAX,bool useEntrez=false, std::string granularity="", std::unordered_map<std::string,std::vector<std::string>> typeToNodeNames = std::unordered_map<std::string,std::vector<std::string>>(), bool undirectedTypeEdges = false);
/**
 * \brief   Returns the new virtual nodes associated with a type along the edges in the augmented graph, it also return the graph of the interactions between types(as a tuple of:
 *         - the start type/agent
 *         - the end type/agent
 *         - the contact times of the interaction, as an unordered set of doubles
 *         granularity needs to be specified as an argument
 * \return the pair (map of the new virtual nodes associated with a type, graph of the interactions between types)
 */
std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>>> interactionContinuousContactsFileToEdgesListAndNodesByName(std::string filename, std::vector<std::string> subtypes,int maximumIntertypeTime=INT32_MAX,bool useEntrez=false, std::string granularity="", std::unordered_map<std::string,std::vector<std::string>> typeToNodeNames = std::unordered_map<std::string,std::vector<std::string>>(), bool undirectedTypeEdges = false, double timestep=1.0);
/**
 * \brief   Returns the new virtual nodes associated with a type along the edges in the augmented graph, uses doubles in principle for contacts it also return the graph of the interactions between types(as a tuple of:
 *         - the start type/agent
 *         - the end type/agent
 *         - the contact times of the interaction, as an unordered set of doubles
 *         granularity needs to be specified as an argument
 * \return the pair (map of the new virtual nodes associated with a type, graph of the interactions between types)
 */
std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>>> interactionContinuousContactsFileToEdgesListAndNodesByName(std::string filename, std::vector<std::string> subtypes,double maximumIntertypeTime=DBL_MAX,bool useEntrez=false, std::string granularity="", std::unordered_map<std::string,std::vector<std::string>> typeToNodeNames = std::unordered_map<std::string,std::vector<std::string>>(), bool undirectedTypeEdges = false, double timestep=1.0);
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
std::map<std::string, std::vector<std::string>> getFullNodesDescription(std::string filename = "resources/graphs/metapathwayNew/nodes.tsv"); 
/**
 * \brief   Return the filenames of all files that have the specified extension
 *          in the specified directory and all subdirectories.
 */
std::vector<std::string> get_all(std::string const & root, std::string const & ext);
/**
 * \brief   Return map of the vector values from the first vector to the second vector
 */
template<typename T>
std::vector<int> get_indexmap_vector_values(std::vector<T> const & origin, std::vector<T> const & toMap){
    std::vector<int> retVec;
    for (int i = 0; i < toMap.size(); ++i) {
        auto it = std::find(origin.begin(), origin.end(), toMap[i]);
        if (it != origin.end()) {
            retVec.push_back(std::distance(origin.begin(), it));
        }
        else {
            std::cout << "[ERROR] utilities::get_indexmap_vector_values : " << toMap[i] << " not found in the origin vector" << std::endl;
            throw std::invalid_argument( "utilities::get_indexmap_vector_values : " + toMap[i] + " not found in the origin vector" );
        }
    }
    return retVec;
}

/**
 * \brief   Return map of the vector values from the second vector to the first vector, if the value is not found in the first vector the value -1 is added to the map
 */
template<typename T>
std::vector<int> get_indexmap_vector_values_full(std::vector<T> const & origin, std::vector<T> const & toMap){
    std::vector<int> retVec;
    std::vector<T> notFoundValuesVector;
    uint notFoundValues = 0;
    for (uint i = 0; i < origin.size(); ++i) {
        auto itorigin = std::find(toMap.begin(), toMap.end(), origin[i]);
        if (itorigin != toMap.end()) {
            retVec.push_back(std::distance(toMap.begin(), itorigin));
        }
        else {
            retVec.push_back(-1);
            notFoundValuesVector.push_back(origin[i]);
            notFoundValues++;
        }
        
    }
    if(notFoundValues == origin.size()){
        std::cout << "[ERROR] utilities::get_indexmap_vector_values_full : all values not found in the origin vector" << std::endl;
        throw std::invalid_argument( "all values not found in the origin vector" );
    }
    if(notFoundValues>0){
        std::cout << "[WARNING] utilities::get_indexmap_vector_values_full : " << notFoundValues << " values not found in the origin vector" << std::endl;
        std::cout << "[WARNING] utilities::get_indexmap_vector_values_full : " << "values not found in the origin vector: ";
        printVector(notFoundValuesVector);
    }
    return retVec;
}

void saveNodeValues(std::string folderName,int iteration, std::string cellName, std::vector<double> nodeValues,std::vector<std::string> nodeNames, bool useEntrez=false, std::string nodesDescriptionFile="");

/**
 * \brief   save node values in folder
 *         better version of the above function
*/
void saveNodeValues(std::string folderName,int iterationOuter, int intraIteration, std::string cellName, std::vector<double> nodeValues,std::vector<std::string> nodeNames, bool useEntrez=false, std::string nodesDescriptionFile="");

/**
 * \brief   save node values in folder
 *         add times as an additional feature
*/
void saveNodeValuesWithTime(std::string folderName,int iterationOuter, int intraIteration, std::string cellName, std::vector<double> nodeValues,std::vector<std::string> nodeNames, bool useEntrez, std::string nodesDescriptionFile="", double timestep=1.0);

/**
 * \brief   save node values in folder, no info about intra-iteration and inter-iteration is passed
 *         add times as an additional feature
*/
void saveNodeValuesWithTimeSimple(std::string folderName, int currentIteration, double currentTime, std::string typeName, std::vector<double> nodeValues,std::vector<std::string> nodeNames, bool useEntrez, std::string nodesDescriptionFile="");


template<typename T>
std::vector<T> vectorSubtraction(std::vector<T> vec1, std::vector<T> vec2){
    if(vec1.size()!=vec2.size()){
        std::cout << "[ERROR] utilities::vectorSubtraction : vectors have different sizes" << std::endl;
        throw std::invalid_argument( "utilities::vectorSubtraction : vectors have different sizes" );
    }
    std::vector<T> retVec;
    for (int i = 0; i < vec1.size(); ++i) {
        retVec.push_back(vec1[i]-vec2[i]);
    }
    return retVec;
}

template<typename T>
std::vector<T> vectorAddition(std::vector<T> vec1, std::vector<T> vec2){
    if(vec1.size()!=vec2.size()){
        std::cout << "[ERROR] utilities::vectorAddition : vectors have different sizes" << std::endl;
        throw std::invalid_argument( "utilities::vectorAddition : vectors have different sizes" );
    }
    std::vector<T> retVec;
    for (int i = 0; i < vec1.size(); ++i) {
        retVec.push_back(vec1[i]+vec2[i]);
    }
    return retVec;
}

template<typename T>
std::vector<T> vectorsIntersection(std::vector<T> vec1, std::vector<T> vec2){
    std::vector<T> retVec;
    for (uint i = 0; i < vec1.size(); ++i) {
        auto it = std::find(vec2.begin(), vec2.end(), vec1[i]);
        if (it != vec2.end()) {
            retVec.push_back(vec1[i]);
        }
    }
    return retVec;
}


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
bool vectorContains(std::vector<T> vec, const T element){
    for (uint i = 0; i < vec.size(); ++i) {
        if(vec[i]==element){
            return true;
        }
    }
    return false;
}