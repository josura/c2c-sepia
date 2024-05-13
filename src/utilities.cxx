#include "utilities.h"
#include <algorithm>
#include <boost/token_functions.hpp>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iterator>
#include <map>
#include <math.h>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;


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

void setRandom(int& val) { 
    int INTMAX = std::numeric_limits<int>::max();
    val =  randomNumber(-INTMAX, INTMAX);
}
void setRandom(double& val) { 
    double DOUBLEMAX = std::numeric_limits<double>::max();
    val = randomRealNumber(-DOUBLEMAX, DOUBLEMAX);
}
void setRandom(char& val) { 
    val = generateRandomCharacter();
}

int SizeToInt(size_t u)
{
    if (u > std::numeric_limits<int>::max())
    {
        throw std::overflow_error(
            "size_t value cannot be stored in a variable of type int.");
    }

    return static_cast<int>(u);
}

long int szudzik(int x, int y)
{
    return x >= y ? x * x + x + y : x + y * y;
}

bool approximatelyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool essentiallyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyGreaterThan(float a, float b, float epsilon)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyLessThan(float a, float b, float epsilon)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool approximatelyEqual(double a, double b, double epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool essentiallyEqual(double a, double b, double epsilon)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyGreaterThan(double a, double b, double epsilon)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyLessThan(double a, double b, double epsilon)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


double hyperbolicTangentScaled(double xInput, double scaleFactor ){
    double firstTerm = std::exp(xInput/scaleFactor);
    double secondTerm = std::exp(-xInput/scaleFactor);
    return scaleFactor*(firstTerm - secondTerm)/(firstTerm + secondTerm);
}

double lerping(double a, double b, double t){
    return a + t * (b - a);
}

bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

bool fileExistsPath(const std::string& filePath)
{
    struct stat info;
    
    // Call the stat function to check if the file exists
    if (stat(filePath.c_str(), &info) != 0) {
        // The stat function returned an error, so the file does not exist
        return false;
    }
    
    // Check if the path corresponds to a regular file
    return (info.st_mode & S_IFREG) != 0;
}

bool folderExists(const std::string& folderPath)
{
    // Check if the folder exists
    bool exists = std::filesystem::is_directory(folderPath);
    
    // Return the result
    return exists;
}

bool createFolder(const std::string& folderPath)
{
    // Check if the folder already exists
    if (folderExists(folderPath)) {
        // The folder already exists
        return true;
    }
    
    // Create the folder
    bool success = std::filesystem::create_directory(folderPath);
    
    // Return the result
    return success;
}

bool setDoubleContainsInterval(std::set<double> set, double lower, double upper){
    if(lower > upper){
        throw std::invalid_argument("utilities::setDoubleContainsInterval: lower bound is greater than upper bound");
    }
    if(set.lower_bound(lower) != set.end() && *set.lower_bound(lower) < upper){
        return true;
    }
    return false;
}

int setDoubleIntervalWidth(std::set<double> set, double lower, double upper){
    if(lower > upper){
        throw std::invalid_argument("utilities::setDoubleContainsInterval: lower bound is greater than upper bound");
    }
    int count = 0;
    for(auto iter = set.cbegin(); iter != set.cend(); iter++){
        if(*iter >= lower && *iter<upper){
            count++;
        }
    }
    return count;
}



std::vector<std::string> splitStringIntoVector(std::string toSplit , std::string delimiter){

    vector<string> tokens;
    boost::split(tokens, toSplit, boost::is_any_of(delimiter));
    return tokens;
}


std::vector<std::string> splitStringIntoVectorTwoParts(std::string toSplit , std::string delimiter){
    std::vector<std::string> ret;
    // split the string into two parts, the first delimiter found is the split point, the delimiter is not included in the result
    int index = toSplit.find(delimiter);
    std::string firstPart = toSplit.substr(0, index);
    std::string secondPart = "";
    if(index >= 0){
        secondPart = toSplit.substr(index + delimiter.length());
    }
    ret.push_back(firstPart);
    if(secondPart.length()>0)
        ret.push_back(secondPart);
    return ret;
}


std::vector<std::string> splitVirtualNodeStringIntoVector(std::string toSplit){
    std::vector<std::string> ret;
    std::vector<std::string> splitted = splitStringIntoVectorTwoParts(toSplit, ":");
    if(splitted.size()==2){
        // first element is the virtual node type (v-in or v-out)
        // second element is the virtual node name (the type of the node, and additionally the name of the node)
        // TODO check if it is a virtual input or a virtual output
        ret.push_back(splitted[0]);
        std::vector<std::string> splittedSecond = splitStringIntoVectorTwoParts(splitted[1], "_");
        if(splittedSecond.size()==2){
            ret.push_back(splittedSecond[0]);
            ret.push_back(splittedSecond[1]);
        } else if(splittedSecond.size()==1){
            ret.push_back(splittedSecond[0]);
        } else {
            throw std::invalid_argument("utilities::splitVirtualNodeStringIntoVector: invalid virtual node string " + toSplit);
        }
    } else {
        throw std::invalid_argument("utilities::splitVirtualNodeStringIntoVector: invalid virtual node string " + toSplit);
    }
    return ret;

}

boost::tokenizer<boost::char_delimiters_separator<char>> splitStringIntoVectorTokenizer(std::string toSplit , char delimiter){
    boost::char_delimiters_separator<char> sep(delimiter);
    boost::tokenizer<boost::char_delimiters_separator<char> > tokens(toSplit, sep);
    return tokens;
}

std::pair<std::vector<int>,std::vector<std::tuple<int,int,double>>> edgesFileToEdgesListByIndex(std::string filename){
    string line;
    std::vector<std::tuple<int,int,double>> ret;
    std::vector<int> nameRet;
    std::unordered_set<int> presentNames;
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==3){
                    int node1 = std::stoi( entries[0]);
                    int node2 = std::stoi(entries[1]);
                    double weight = std::stod( entries[2]);
                    std::tuple<int,int,double> edge(node1,node2,weight);
                    ret.push_back(edge);
                    if(!presentNames.contains(node1)){
                        nameRet.push_back(node1);
                        presentNames.insert(node1);
                    }
                    if(!presentNames.contains(node2)){
                        nameRet.push_back(node2);
                        presentNames.insert(node2);

                    }
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::edgeFileEdgesListByIndex: file does not exists");
    }
    return std::pair<std::vector<int>,std::vector<std::tuple<int,int,double>>> (nameRet,ret);
}

std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>> edgesFileToEdgesListAndNodesByName(std::string filename){
    string line;
    std::vector<std::tuple<std::string,std::string,double>> ret;
    std::vector<std::string> nameRet;
    std::unordered_set<std::string> presentNames;
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> entriesHeader = splitStringIntoVector(line, "\t");
            int indexStart=-1,indexEnd=-1,indexWeight=-1;
            for(uint i = 0; i < entriesHeader.size(); i++){
                if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("start") != std::string::npos) {
                    indexStart = i;
                }
                else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("end") != std::string::npos) {
                    indexEnd = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("weight") != std::string::npos) {
                    indexWeight = i;
                }
            }

            if(indexStart < 0 || indexEnd < 0 || indexWeight < 0){
                if(entriesHeader.size()==3){
                    indexStart = 0;
                    indexEnd = 1;
                    indexWeight = 2;
                    std::cout << "[WARNING] using the first, second and third column as start, end and weight in the graph file:" << filename << std::endl;
                } else {
                    std::string error = "utilities::edgesFileToEdgesListAndNodesByName: header of file" + filename + " does not contain start, end or weight";
                    throw std::invalid_argument(error);
                }
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==entriesHeader.size()){
                    std::string node1 = entries[indexStart];
                    std::string node2 = entries[indexEnd];
                    double weight = std::stod( entries[indexWeight]);
                    std::tuple<std::string,std::string,double> edge(node1,node2,weight);
                    ret.push_back(edge);
                    if(!presentNames.contains(node1)){
                        nameRet.push_back(node1);
                        presentNames.insert(node1);
                    }
                    if(!presentNames.contains(node2)){
                        nameRet.push_back(node2);
                        presentNames.insert(node2);

                    }
                }
            }
            // control if resulting edges vector is empty
            if(ret.size()==0){
                std::cerr << "[WARNING] edgesFileToEdgesListAndNodesByName: no edges found in the file " << filename << " .Use the nodeDescriptionFolder parameter to pass the graphs nodes, otherwise an error will occur" << std::endl;
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::edgeFileEdgesListByIndex: file does not exists " + filename);
    }
    return std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>> (nameRet,ret);
}

std::pair<std::vector<std::string>,std::vector<std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>>>> edgesFileToEdgesListAndNodesByNameFromFolder(std::string filename){
    std::vector<std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>>> ret;
    std::vector<std::string> graphNames;
    std::vector<std::string> files = get_all(filename,".tsv");
    for(auto iter = files.cbegin(); iter!=files.cend();iter++){
        std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>> edgesListAndNodesByName = edgesFileToEdgesListAndNodesByName(*iter);
        std::vector<std::string> splitted = splitStringIntoVector(*iter, "/"); //split the path
        std::string filename = splitted[splitted.size()-1]; //last element
        std::vector<std::string> splittedFilename = splitStringIntoVector(filename, "."); //split the extension
        graphNames.push_back(splittedFilename[0]);
        ret.push_back(edgesListAndNodesByName);
    }
    return std::pair<std::vector<std::string>,std::vector<std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>>>>(graphNames,ret);
}

std::vector<std::string> getTypesFromFolderFileNames(std::string folderPath){
    std::vector<std::string> ret;
    std::vector<std::string> files = get_all(folderPath,".tsv");
    for(auto iter = files.cbegin(); iter!=files.cend();iter++){
        std::vector<std::string> splitted = splitStringIntoVector(*iter, "/"); //split the path
        std::string filename = splitted[splitted.size()-1]; //last element
        std::vector<std::string> splittedFilename = splitStringIntoVector(filename, "."); //split the extension
        ret.push_back(splittedFilename[0]);
    }
    return ret;
}

std::vector<std::string> getTypesFromMatrixFile(std::string matrixFilepath){
    std::vector<std::string> typeNames;
    std::string line;
    if(file_exists(matrixFilepath)){
        ifstream myfile (matrixFilepath);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> splittedHeader = splitStringIntoVector(line, "\t");  //could already be used as the cellnames vector,
            for (int i = 1; i < SizeToInt( splittedHeader.size()); i++) {
                typeNames.push_back(splittedHeader[i]);
            }
            myfile.close();
        }
        //check for duplicates
        if(controlForDuplicates(typeNames)){
            throw std::invalid_argument("utilities::getTypesFromMatrixFile: duplicate types in the matrix file, aborting " + matrixFilepath);
        }
    } else {
        throw std::invalid_argument("utilities::getTypesFromMatrixFile: file does not exists " + matrixFilepath);
    }
    return typeNames;
}

std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeMatrixToCellVectors(std::string filename, const std::vector<std::string>& finalNames,bool useEntrez){
    string line;
    std::vector<std::vector<double>> ret;
    std::vector<std::string> cellNames;
    std::vector<std::string> geneNames;
    std::vector<std::string> discardedGenes;
    std::map<std::string, int> finalGenesToIndex;
    for(int i = 0 ; i < SizeToInt(finalNames.size()); i++){
        finalGenesToIndex[finalNames[i]] = i;
    }
    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> splittedHeader = splitStringIntoVector(line, "\t");  //could already be used as the cellnames vector,
            for (int i = 1; i < SizeToInt( splittedHeader.size()); i++) {
                cellNames.push_back(splittedHeader[i]);
                ret.push_back(std::vector<double>(finalNames.size(),0));
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==splittedHeader.size()){
                    if(!useEntrez){
                        geneNames.push_back(entries[0]);
                        for(int i = 1; i < SizeToInt(entries.size());i++){
                            //ret[i-1].push_back(std::stod(entries[i]));
                            ret[i-1][finalGenesToIndex[entries[0]]] = std::stod(entries[i]);
                        } //TODO control over the genes in the metapathway like its done below with the mapping, but without the mapping and by taking a vector maybe
                    }
                    else{
                        if (mapEnsembleToEntrez.contains(entries[0]) && finalGenesToIndex.contains(mapEnsembleToEntrez[entries[0]])) {
                            geneNames.push_back(mapEnsembleToEntrez[entries[0]]);
                            for(int i = 1; i < SizeToInt(entries.size());i++){
                                //ret[i-1].push_back(std::stod(entries[i]));
                                ret[i-1][finalGenesToIndex[mapEnsembleToEntrez[entries[0]]]] = std::stod(entries[i]);
                            }
                        } else{
                            discardedGenes.push_back(entries[0]);
                        }//else don't do nothing since the node is not in the graph
                    }
                } else {
                    throw std::invalid_argument("utilities::edgeFileEdgesListByIndex: header doesn't have the same amount of columns as the data " + filename);
                }
            }
            myfile.close();
            std::cout << "[LOG] No node in the type specified for nodes: " << std::endl;
            for(auto iter = discardedGenes.cbegin();iter!=discardedGenes.cend();iter++){
                std::cout << "," << *iter;
            }
            std::cout << std::endl <<"[LOG] discarding values for the nodes not in the graph" << std::endl;
            
        }
    } else {
        throw std::invalid_argument("utilities::edgeFileEdgesListByIndex: file does not exists " + filename);
    }
    return std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> (geneNames,cellNames,ret);

}


// TODO: refactor
std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeMatrixToCellVectors(std::string filename, const std::vector<std::string>& finalNames,std::vector<std::string> subTypes ,bool useEntrez){
    string line;
    std::vector<std::vector<double>> ret;
    std::vector<std::string> cellNames;
    std::vector<std::string> geneNames;
    std::vector<std::string> discardedGenes;
    std::map<std::string, int> finalGenesToIndex;
    for(int i = 0 ; i < SizeToInt(finalNames.size()); i++){
        finalGenesToIndex[finalNames[i]] = i;
    }
    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> splittedHeader = splitStringIntoVector(line, "\t");  //could already be used as the cellnames vector,
            std::vector<int> subTypeIndexes;
            for (int i = 1; i < SizeToInt( splittedHeader.size()); i++) {
                if(vectorContains(subTypes,splittedHeader[i])){
                    cellNames.push_back(splittedHeader[i]);
                    subTypeIndexes.push_back(i);
                    ret.push_back(std::vector<double>(finalNames.size(),0));
                }
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==splittedHeader.size()){
                    if(!useEntrez){
                        geneNames.push_back(entries[0]);
                        for(uint i = 0; i<subTypeIndexes.size();i++){
                            ret[i][finalGenesToIndex[entries[0]]] = std::stod(entries[subTypeIndexes[i]]);
                        } //TODO control over the genes in the metapathway like its done below with the mapping, but without the mapping and by taking a vector maybe
                    }
                    else{
                        if (mapEnsembleToEntrez.contains(entries[0]) && finalGenesToIndex.contains(mapEnsembleToEntrez[entries[0]])) {
                            geneNames.push_back(mapEnsembleToEntrez[entries[0]]);
                            for(uint i = 0; i< subTypeIndexes.size();i++){
                                ret[i][finalGenesToIndex[mapEnsembleToEntrez[entries[0]]]] = std::stod(entries[subTypeIndexes[i]]);
                            }
                        } else{
                            discardedGenes.push_back(entries[0]);
                        }//else don't do nothing since the node is not in the graph
                    }
                } else {
                    throw std::invalid_argument("utilities::edgeFileEdgesListByIndex: header doesn't have the same amount of columns as the data " + filename);
                }
            }
            myfile.close();
            std::cout << "[LOG] No nodes in the graph for nodes: " << std::endl;
            for(auto iter = discardedGenes.cbegin();iter!=discardedGenes.cend();iter++){
                std::cout << "," << *iter;
            }
            std::cout << std::endl <<"[LOG] discarding values for the nodes not in the graph" << std::endl;
            
        }
    } else {
        throw std::invalid_argument("utilities::edgeFileEdgesListByIndex: file does not exists " + filename);
    }
    return std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> (geneNames,cellNames,ret);

}

// TODO: refactor
std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeCellVectorsFromFolder(std::string folderPath,const std::vector<std::string>& allTypes , const std::vector<std::vector<std::string>>& finalNames,std::vector<std::string> subType, bool useEntrez){
    std::vector<std::string> cellNames;
    std::vector<std::string> geneNames;
    std::vector<std::vector<double>> ret;
    std::vector<std::string> discardedGenes;
    std::map<std::string ,std::map<std::string, int>> finalGenesToIndex;
    for(int i = 0 ; i < SizeToInt(finalNames.size()); i++){
        std::map<std::string, int> tmp;
        for(int j = 0 ; j < SizeToInt(finalNames[i].size()); j++){
            tmp[finalNames[i][j]] = j;
        }
        finalGenesToIndex[allTypes[i]]= tmp;
    }
    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    auto files = get_all(folderPath,".tsv");
    if(files.size()==0){
        throw std::invalid_argument("utilities::logFoldChangeCellVectorsFromFolder: no files found in the folder " + folderPath);
    }
    //default argument for subtype is empty, if empty, use all files in the folder
    if(subType.size()==0){
        for(auto iter = files.cbegin();iter!=files.cend();iter++){
            std::vector<std::string> splitted = splitStringIntoVector(*iter, "/"); //split the path
            std::string filename = splitted[splitted.size()-1]; //last element
            subType.push_back(splitStringIntoVector(filename, ".")[0]);
        }
    }
    //filter files from subtypes (first part of the filename before the extension)
    std::vector<std::string> filteredFiles;
    for(auto iter = files.cbegin();iter!=files.cend();iter++){
        std::vector<std::string> fileSplitPath = splitStringIntoVector(*iter, "/");
        std::string filename = fileSplitPath[fileSplitPath.size()-1];
        std::string type = splitStringIntoVector(filename, ".")[0];
        if(vectorContains(subType,type)){
            filteredFiles.push_back(*iter);
        } else {
            std::cout << "[LOG] discarding file " << *iter << " since it is not in the subtypes" << std::endl;
        }
    }
    if(filteredFiles.size()==0){
        throw std::invalid_argument("utilities::logFoldChangeCellVectorsFromFolder: no files found in the folder that are similar to the subtypes " + folderPath);
    }
    for(auto iter = filteredFiles.cbegin();iter!=filteredFiles.cend();iter++){
        std::vector<std::string> splitted = splitStringIntoVector(*iter, "/"); //split the path
        std::string filenameNoPath = splitted[splitted.size()-1]; //last element
        std::string cellName = splitStringIntoVector(filenameNoPath, ".")[0];
        cellNames.push_back(cellName);
        std::string filename = *iter;
        if(file_exists(filename)){
            //first line is the header, the first column is the gene, the second column is the value
            ifstream myfile (filename);
            string line;
            std::vector<double> cellValues(finalGenesToIndex[cellName].size(),0);
            std::string lineHeader;
            getline (myfile,lineHeader);  // first line is header IMPORTANT
            std::vector<std::string> splittedHeader = splitStringIntoVector(lineHeader, "\t");
            //check if the header is correct
            if(splittedHeader.size()!=2){
                throw std::invalid_argument("utilities::logFoldChangeCellVectorsFromFolder: header doesn't have the same amount of columns as the data " + filename);
            }
            if(splittedHeader[0]!="name" || splittedHeader[1] != "value"){
                throw std::invalid_argument("utilities::logFoldChangeCellVectorsFromFolder: header doesn't have the name and value columns or it does not have  an header" + filename);
            }
            //get file contents
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==2){
                    if(!useEntrez){
                        if(finalGenesToIndex[cellName].contains(entries[0])){
                            cellValues[finalGenesToIndex[cellName][entries[0]]] = std::stod(entries[1]);
                            geneNames.push_back(entries[0]);
                        } else{
                            discardedGenes.push_back(entries[0]);
                        }
                    }
                    else{
                        if (mapEnsembleToEntrez.contains(entries[0]) && finalGenesToIndex[cellName].contains(mapEnsembleToEntrez[entries[0]])) {
                            cellValues[finalGenesToIndex[cellName][mapEnsembleToEntrez[entries[0]]]] = std::stod(entries[1]);
                            geneNames.push_back(mapEnsembleToEntrez[entries[0]]);
                        } else{
                            discardedGenes.push_back(entries[0]);
                        }//else don't do nothing since the node is not in the graph
                    }
                } else {
                    throw std::invalid_argument("utilities::logFoldChangeCellVectorsFromFolder: header doesn't have the same amount of columns as the data " + filename);
                }
            }
            myfile.close();
            std::cout << "[LOG] discarding values for the nodes not in the graph for type "<< cellName << ", the nodes discarded are:" << std::endl;
            for(auto iter = discardedGenes.cbegin();iter!=discardedGenes.cend();iter++){
                std::cout << "," << *iter;
            }
            std::cout << std::endl;
            ret.push_back(cellValues);

        
        } else {
            throw std::invalid_argument("utilities::logFoldChangeCellVectorsFromFolder: file does not exists " + filename);
        }
    }
    return std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> (geneNames,cellNames,ret);
            
}


std::map<std::string,std::vector<std::string>> nodeNamesFromFolder(std::string folderPath){
    std::map<std::string,std::vector<std::string>> ret;
    auto files = get_all(folderPath,".tsv");
    if(files.size()==0){
        throw std::invalid_argument("utilities::nodeNamesFromFolder: no files found in the folder " + folderPath);
    }
    for(auto iter = files.cbegin();iter!=files.cend();iter++){
        std::vector<std::string> splitted = splitStringIntoVector(*iter, "/"); //split the path
        std::string filename = splitted[splitted.size()-1]; //last element
        std::vector<std::string> splittedFilename = splitStringIntoVector(filename, "."); //split the extension
        std::string type = splittedFilename[0];
        std::vector<std::string> nodeNames;
        if(file_exists(*iter)){
            ifstream myfile (*iter);
            string line;
            getline (myfile,line);  // first line is header
            // search for name column
            std::vector<std::string> splittedHeader = splitStringIntoVector(line, "\t");
            int indexName=-1;
            for(uint i = 0; i < splittedHeader.size(); i++){
                if (boost::algorithm::to_lower_copy(splittedHeader[i]).find("name") != std::string::npos) {
                    indexName = i;
                }
            }
            if(indexName < 0){
                throw std::invalid_argument("utilities::nodeNamesFromFolder: invalid file, the header does not contain a name feature");
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==splittedHeader.size()){
                    nodeNames.push_back(entries[indexName]);
                }
            }
            myfile.close();
            ret[type] = nodeNames;
        } else {
            throw std::invalid_argument("utilities::nodeNamesFromFolder: file does not exists " + *iter);
        }
    }
    return ret;

}



std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> interactionFileToEdgesListAndNodesByName(std::string filename,bool useEntrez){
    string line;
    std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> ret;
    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> entriesHeader = splitStringIntoVector(line, "\t");
            int indexTypeStart=-1,indexTypeEnd=-1,indexStartNode=-1,indexEndNode=-1,indexWeight=-1;
            //TODO change attributes names to be more general
            for(uint i = 0; i < entriesHeader.size(); i++){
                if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("starttype") != std::string::npos) {
                    indexTypeStart = i;
                }
                else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endtype") != std::string::npos) {
                    indexTypeEnd = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("startnodename") != std::string::npos) {
                    indexStartNode = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endnodename") != std::string::npos) {
                    indexEndNode = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("weight") != std::string::npos) {
                    indexWeight = i;
                }
                //startType	startNodeName	endType	endNodeName	weight
            }
            if(indexTypeStart < 0 || indexTypeEnd < 0 || indexStartNode < 0 || indexEndNode < 0 || indexWeight < 0){
                throw std::invalid_argument("utilities::interactionFileToEdgesListAndNodesByName: invalid file, the header does not contain a startType, or an endType, or a startNodeName, or a endNodeName, or a weight feature");
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==5){
                    std::string startNodeName,endNodeName;
                    if(!useEntrez){
                        startNodeName = entries[indexStartNode];
                        endNodeName = entries[indexEndNode];    

                        std::string startType = entries[indexTypeStart];
                        std::string endType = entries[indexTypeEnd];
                        double weight = std::stod( entries[indexWeight]);
                        std::string virtualInputEndType = "v-in:" + startType;
                        std::string virtualOutputstartType = "v-out:" + endType;
                        std::tuple<std::string,std::string,double> edgestartType(startNodeName, virtualOutputstartType,weight);
                        std::tuple<std::string,std::string,double> edgeEndType(virtualInputEndType, endNodeName,weight);
                        if(ret.contains(startType)){
                            ret[startType].push_back(edgestartType);
                        }else{
                            ret[startType] = std::vector<std::tuple<std::string,std::string,double>>();
                            ret[startType].push_back(edgestartType);
                        }

                        if(ret.contains(endType)){
                            ret[endType].push_back(edgeEndType);
                        }else{
                            ret[endType] = std::vector<std::tuple<std::string,std::string,double>>();
                            ret[endType].push_back(edgeEndType);
                        }
                    } else{
                        if(mapEnsembleToEntrez.contains(entries[indexStartNode]) && mapEnsembleToEntrez.contains(entries[indexEndNode])){
                            startNodeName = mapEnsembleToEntrez[entries[indexStartNode]];
                            endNodeName = mapEnsembleToEntrez[entries[indexEndNode]];

                            std::string startType = entries[indexTypeStart];
                            std::string endType = entries[indexTypeEnd];
                            double weight = std::stod( entries[indexWeight]);
                            std::string virtualInputEndType = "v-in:" + startType;
                            std::string virtualOutputstartType = "v-out:" + endType;
                            std::tuple<std::string,std::string,double> edgestartType(startNodeName, virtualOutputstartType,weight);
                            std::tuple<std::string,std::string,double> edgeEndType(virtualInputEndType, endNodeName,weight);
                            if(ret.contains(startType)){
                                ret[startType].push_back(edgestartType);
                            }else{
                                ret[startType] = std::vector<std::tuple<std::string,std::string,double>>();
                                ret[startType].push_back(edgestartType);
                            }

                            if(ret.contains(endType)){
                                ret[endType].push_back(edgeEndType);
                            }else{
                                ret[endType] = std::vector<std::tuple<std::string,std::string,double>>();
                                ret[endType].push_back(edgeEndType);
                            }
                        }
                    }
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::interactionFileToEdgesListAndNodesByName: file does not exists " + filename);
    }
    return ret;
}

std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> interactionFileToEdgesListAndNodesByName(std::string filename,std::vector<std::string> subtypes,bool useEntrez){
    string line;
    std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> ret;
    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> entriesHeader = splitStringIntoVector(line, "\t");
            int indexTypeStart=-1,indexTypeEnd=-1,indexStartNode=-1,indexEndNode=-1,indexWeight=-1;
            for(uint i = 0; i < entriesHeader.size(); i++){ //TODO change names in the header to be more general
                if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("starttype") != std::string::npos) {
                    indexTypeStart = i;
                }
                else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endtype") != std::string::npos) {
                    indexTypeEnd = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("startnodename") != std::string::npos) {
                    indexStartNode = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endnodename") != std::string::npos) {
                    indexEndNode = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("weight") != std::string::npos) {
                    indexWeight = i;
                }
                //startType	startNodeName	endType	endNodeName	weight
            }
            if(indexTypeStart < 0 || indexTypeEnd < 0 || indexStartNode < 0 || indexEndNode < 0 || indexWeight < 0){
                throw std::invalid_argument("utilities::interactionFileToEdgesListAndNodesByName: invalid file, the header does not contain a startType, or an endType, or a Ligand gene, or a receptor gene, or a weight feature");
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==5){
                    std::string startNodeName,endNodeName;
                    if(!useEntrez){
                        startNodeName = entries[indexStartNode];
                        endNodeName = entries[indexEndNode];    

                        std::string startType = entries[indexTypeStart];
                        std::string endType = entries[indexTypeEnd];
                        if(vectorContains(subtypes, startType) && vectorContains(subtypes, endType)){
                            double weight = std::stod( entries[indexWeight]);
                            std::string virtualInputEndType = "v-in:" + startType;
                            std::string virtualOutputstartType = "v-out:" + endType;
                            std::tuple<std::string,std::string,double> edgestartType(startNodeName, virtualOutputstartType,weight);
                            std::tuple<std::string,std::string,double> edgeEndType(virtualInputEndType, endNodeName,weight);
                            if(ret.contains(startType)){
                                ret[startType].push_back(edgestartType);
                            }else{
                                ret[startType] = std::vector<std::tuple<std::string,std::string,double>>();
                                ret[startType].push_back(edgestartType);
                            }

                            if(ret.contains(endType)){
                                ret[endType].push_back(edgeEndType);
                            }else{
                                ret[endType] = std::vector<std::tuple<std::string,std::string,double>>();
                                ret[endType].push_back(edgeEndType);
                            }
                        } else {
                            //ignored because not in the subtypes
                        }
                        
                    } else{
                        if(mapEnsembleToEntrez.contains(entries[indexStartNode]) && mapEnsembleToEntrez.contains(entries[indexEndNode])){
                            startNodeName = mapEnsembleToEntrez[entries[indexStartNode]];
                            endNodeName = mapEnsembleToEntrez[entries[indexEndNode]];

                            std::string startType = entries[indexTypeStart];
                            std::string endType = entries[indexTypeEnd];
                            if(vectorContains(subtypes, startType) && vectorContains(subtypes, endType)){
                                double weight = std::stod( entries[indexWeight]);
                                std::string virtualInputEndType = "v-in:" + startType;
                                std::string virtualOutputstartType = "v-out:" + endType;
                                std::tuple<std::string,std::string,double> edgestartType(startNodeName, virtualOutputstartType,weight);
                                std::tuple<std::string,std::string,double> edgeEndType(virtualInputEndType, endNodeName,weight);
                                if(ret.contains(startType)){
                                    ret[startType].push_back(edgestartType);
                                }else{
                                    ret[startType] = std::vector<std::tuple<std::string,std::string,double>>();
                                    ret[startType].push_back(edgestartType);
                                }

                                if(ret.contains(endType)){
                                    ret[endType].push_back(edgeEndType);
                                }else{
                                    ret[endType] = std::vector<std::tuple<std::string,std::string,double>>();
                                    ret[endType].push_back(edgeEndType);
                                }
                            } else {
                                //ignored because not in the subtypes
                            }
                        }
                    }
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::interactionFileToEdgesListAndNodesByName: file does not exists " + filename);
    }
    return ret;
}

std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::unordered_set<int>, double>>> interactionContactsFileToEdgesListAndNodesByName(std::string filename, std::vector<std::string> subtypes, int maximumIntertypeTime, bool useEntrez, std::string granularity,std::unordered_map<std::string,std::vector<std::string>> typeToNodeNames, bool undirectedTypeEdges){
    string line;
    std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::unordered_set<int>, double>>> ret;
    // control if the granularity is valid
    if(granularity != "" && granularity != "type" && granularity != "node" && granularity != "typeAndNode"){
        throw std::invalid_argument("utilities::interactionContactsFileToEdgesListAndNodesByName: invalid granularity, it must be typeAndNode(finer) or type(coarser), or only node(no types)");
    }
    if(granularity == ""){
        granularity = "type";
    }
    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> entriesHeader = splitStringIntoVector(line, "\t");
            int indexTypeStart=-1, indexTypeEnd=-1, indexStartNode=-1, indexEndNode=-1, indexWeight=-1, indexContactTimes=-1;
            for(uint i = 0; i < entriesHeader.size(); i++){ //TODO change names in the header to be more general
                if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("starttype") != std::string::npos) {
                    indexTypeStart = i;
                }
                else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endtype") != std::string::npos) {
                    indexTypeEnd = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("startnodename") != std::string::npos) {
                    indexStartNode = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endnodename") != std::string::npos) {
                    indexEndNode = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("weight") != std::string::npos) {
                    indexWeight = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("contacttimes") != std::string::npos) {
                    indexContactTimes = i;
                }
                //startType	startNodeName	endType	endNodeName	weight
            }
            bool noContactTimes = false;
            if(indexTypeStart < 0 || indexTypeEnd < 0 || indexStartNode < 0 || indexEndNode < 0 || indexWeight < 0){
                throw std::invalid_argument("utilities::interactionFileToEdgesListAndNodesByName: invalid file, the header does not contain a startType, or an endType, or a start node, or an end node, or a weight feature");
            }
            if(indexContactTimes < 0){
                noContactTimes = true;
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==5 || entries.size()==6){
                    std::string startNodeName,endNodeName;
                    if(!useEntrez){
                        startNodeName = entries[indexStartNode];
                        endNodeName = entries[indexEndNode];
                        
                    } else{
                        if(mapEnsembleToEntrez.contains(entries[indexStartNode]) && mapEnsembleToEntrez.contains(entries[indexEndNode])){
                            startNodeName = mapEnsembleToEntrez[entries[indexStartNode]];
                            endNodeName = mapEnsembleToEntrez[entries[indexEndNode]];
                        }
                    }
                    std::string startType = entries[indexTypeStart];
                    std::string endType = entries[indexTypeEnd];
                    if(typeToNodeNames.size() != 0){
                        if(!vectorContains(typeToNodeNames[startType],startNodeName)){
                            std::cout << "[ERROR] start node <"<< startNodeName <<"> for type: " << startType << " is not in the specified network, aborting " <<std::endl;
                            throw std::invalid_argument("utilities::interactionContactsFileToEdgesListAndNodesByName: invalid file, the start node" + startNodeName + " is not in the type specified, aborting");
                        }

                        if(!vectorContains(typeToNodeNames[endType],endNodeName)){
                            std::cout << "[ERROR] end node <"<< endNodeName <<"> for type: " << endType << " is not in the specified network, aborting " <<std::endl;
                            throw std::invalid_argument("utilities::interactionContactsFileToEdgesListAndNodesByName: invalid file, the end node" + endNodeName + " is not in the type specified, aborting");
                        }
                    }

                    std::unordered_set<int> contactTimes;
                    double weight = std::stod( entries[indexWeight]);
                    if(noContactTimes){
                        // if no contact times are specified, then every time is a contact time, so all contact times from 0 to maximumIntertypeTime are added
                        for(int i = 0; i < maximumIntertypeTime; i++){
                            contactTimes.insert(i);
                        }
                    } else {
                        std::string contactTimesString = entries[indexContactTimes];
                        std::vector<std::string> splittedContactTimes = splitStringIntoVector(contactTimesString, ",");
                        for(auto iter = splittedContactTimes.cbegin(); iter!=splittedContactTimes.cend(); iter++){
                            int contactTime = std::stoi(*iter);
                            if(contactTime <= maximumIntertypeTime){
                                contactTimes.insert(contactTime);
                            }
                        }
                    }
                    // add the edge to the augmented graph, only if the startType and the endType are in the subtypes
                    if(vectorContains(subtypes, startType) && vectorContains(subtypes, endType)){
                        std::string virtualInputEndType = "";
                        std::string virtualOutputStartType = "";
                        // when undirectedTypeEdges is true, virtual input for the startType and virtual output for the endType are added to the respective graphs 
                        std::string virtualInputStartType = "";
                        std::string virtualOutputEndType = "";
                        if(granularity == "typeAndNode"){
                            virtualInputEndType = "v-in:" + startType + "_" + startNodeName;
                            virtualOutputStartType = "v-out:" + endType + "_" + endNodeName;
                            virtualInputStartType = "v-in:" + endType + "_" + endNodeName;
                            virtualOutputEndType = "v-out:" + startType + "_" + startNodeName;
                        } else if (granularity == "type"){
                            virtualInputEndType = "v-in:" + startType;
                            virtualOutputStartType = "v-out:" + endType;
                            virtualInputStartType = "v-in:" + endType;
                            virtualOutputEndType = "v-out:" + startType;
                        } else {
                            virtualInputEndType = "v-in:" + startNodeName;
                            virtualOutputStartType = "v-out:" + endNodeName;
                            virtualInputStartType = "v-in:" + endNodeName;
                            virtualOutputEndType = "v-out:" + startNodeName;
                        }
                        std::tuple<std::string,std::string,double> edgestartType(startNodeName, virtualOutputStartType,weight);
                        std::tuple<std::string,std::string,double> undirectedEdgestartType(virtualInputStartType, startNodeName,weight);
                        std::tuple<std::string,std::string,double> edgeEndType(virtualInputEndType, endNodeName,weight);
                        std::tuple<std::string,std::string,double> undirectedEdgeEndType(endNodeName, virtualOutputEndType,weight);
                        // add the edge to the startType
                        if(!ret.first.contains(startType)){
                            ret.first[startType] = std::vector<std::tuple<std::string,std::string,double>>();
                        }
                        ret.first[startType].push_back(edgestartType);
                        if(undirectedTypeEdges){
                            ret.first[startType].push_back(undirectedEdgestartType);
                        }

                        // add the edge to the endType
                        if(!ret.first.contains(endType)){
                            ret.first[endType] = std::vector<std::tuple<std::string,std::string,double>>();
                        }
                        ret.first[endType].push_back(edgeEndType);
                        if(undirectedTypeEdges){
                            ret.first[endType].push_back(undirectedEdgeEndType);
                        }
                    } else {
                        //ignored because not in the subtypes
                    }
                    // add the edge with the contact times to the second vector in ret
                    ret.second.push_back(std::tuple<std::string, std::string, std::string, std::string, std::unordered_set<int>, double>(startNodeName, endNodeName, startType, endType, contactTimes, weight));
                    if(undirectedTypeEdges){
                        ret.second.push_back(std::tuple<std::string, std::string, std::string, std::string, std::unordered_set<int>, double>(endNodeName, startNodeName, endType, startType, contactTimes, weight));
                    }
                } else {
                    std::cout << "[ERROR] columns detected: " << entries.size() << " columns " <<std::endl;
                    throw std::invalid_argument("utilities::interactionFileToEdgesListAndNodesByName: header doesn't have the right amount of columns(5 or 6 when considering interaction times) ");
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::interactionFileToEdgesListAndNodesByName: file does not exists " + filename);
    }
    return ret;
}

// considering a timestep, the function is almost the same as above but generate double contact times, based upon the timestep
std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>>> interactionContinuousContactsFileToEdgesListAndNodesByName(std::string filename, std::vector<std::string> subtypes, int maximumIntertypeTime, bool useEntrez, std::string granularity,std::unordered_map<std::string,std::vector<std::string>> typeToNodeNames, bool undirectedTypeEdges, double timestep){
    string line;
    std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>>> ret;
    // control if the granularity is valid
    if(granularity != "" && granularity != "type" && granularity != "node" && granularity != "typeAndNode"){
        throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: invalid granularity, it must be typeAndNode(finer) or type(coarser), or only node(no types)");
    }
    if(granularity == ""){
        granularity = "type";
    }

    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> entriesHeader = splitStringIntoVector(line, "\t");
            int indexTypeStart=-1, indexTypeEnd=-1, indexStartNode=-1, indexEndNode=-1, indexWeight=-1, indexContactTimes=-1;
            for(uint i = 0; i < entriesHeader.size(); i++){ //TODO change names in the header to be more general
                if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("starttype") != std::string::npos) {
                    indexTypeStart = i;
                }
                else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endtype") != std::string::npos) {
                    indexTypeEnd = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("startnodename") != std::string::npos) {
                    indexStartNode = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endnodename") != std::string::npos) {
                    indexEndNode = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("weight") != std::string::npos) {
                    indexWeight = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("contacttimes") != std::string::npos) {
                    indexContactTimes = i;
                }
                //startType	startNodeName	endType	endNodeName	weight
            }
            bool noContactTimes = false;
            if(indexTypeStart < 0 || indexTypeEnd < 0 || indexStartNode < 0 || indexEndNode < 0 || indexWeight < 0){
                throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: invalid file, the header does not contain a startType, or an endType, or a start node, or an end node, or a weight feature");
            }
            if(indexContactTimes < 0){
                noContactTimes = true;
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==5 || entries.size()==6){
                    std::string startNodeName,endNodeName;
                    if(!useEntrez){
                        startNodeName = entries[indexStartNode];
                        endNodeName = entries[indexEndNode];
                        
                    } else{
                        if(mapEnsembleToEntrez.contains(entries[indexStartNode]) && mapEnsembleToEntrez.contains(entries[indexEndNode])){
                            startNodeName = mapEnsembleToEntrez[entries[indexStartNode]];
                            endNodeName = mapEnsembleToEntrez[entries[indexEndNode]];
                        }
                    }
                    std::string startType = entries[indexTypeStart];
                    std::string endType = entries[indexTypeEnd];
                    if(typeToNodeNames.size() != 0){
                        if(!vectorContains(typeToNodeNames[startType],startNodeName)){
                            std::cout << "[ERROR] start node <"<< startNodeName <<"> for type: " << startType << " is not in the specified network, aborting " <<std::endl;
                            throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: invalid file, the start node" + startNodeName + " is not in the type specified, aborting");
                        }

                        if(!vectorContains(typeToNodeNames[endType],endNodeName)){
                            std::cout << "[ERROR] end node <"<< endNodeName <<"> for type: " << endType << " is not in the specified network, aborting " <<std::endl;
                            throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: invalid file, the end node" + endNodeName + " is not in the type specified, aborting");
                        }
                    }

                    std::set<double> contactTimes;
                    double weight = std::stod( entries[indexWeight]);
                    if(noContactTimes){
                        // if no contact times are specified, then every time is a contact time, so all contact times from 0 to maximumIntertypeTime are added
                        for(int i = 0; i < maximumIntertypeTime; i++){
                            contactTimes.insert(i*timestep);
                        }
                    } else {
                        std::string contactTimesString = entries[indexContactTimes];
                        std::vector<std::string> splittedContactTimes = splitStringIntoVector(contactTimesString, ",");
                        for(auto iter = splittedContactTimes.cbegin(); iter!=splittedContactTimes.cend(); iter++){
                            double contactTime = std::stod(*iter);
                            // TODO change interTypeTime to be a double since it is a time, and I am now using the timestep to define the times and quantize the contacts
                            if(contactTime <= maximumIntertypeTime){
                                contactTimes.insert(contactTime);
                            }
                        }
                    }
                    // add the edge to the augmented graph, only if the startType and the endType are in the subtypes
                    if(vectorContains(subtypes, startType) && vectorContains(subtypes, endType)){
                        std::string virtualInputEndType = "";
                        std::string virtualOutputStartType = "";
                        // when undirectedTypeEdges is true, virtual input for the startType and virtual output for the endType are added to the respective graphs 
                        std::string virtualInputStartType = "";
                        std::string virtualOutputEndType = "";
                        if(granularity == "typeAndNode"){
                            virtualInputEndType = "v-in:" + startType + "_" + startNodeName;
                            virtualOutputStartType = "v-out:" + endType + "_" + endNodeName;
                            virtualInputStartType = "v-in:" + endType + "_" + endNodeName;
                            virtualOutputEndType = "v-out:" + startType + "_" + startNodeName;
                        } else if (granularity == "type"){
                            virtualInputEndType = "v-in:" + startType;
                            virtualOutputStartType = "v-out:" + endType;
                            virtualInputStartType = "v-in:" + endType;
                            virtualOutputEndType = "v-out:" + startType;
                        } else {
                            virtualInputEndType = "v-in:" + startNodeName;
                            virtualOutputStartType = "v-out:" + endNodeName;
                            virtualInputStartType = "v-in:" + endNodeName;
                            virtualOutputEndType = "v-out:" + startNodeName;
                        }
                        std::tuple<std::string,std::string,double> edgestartType(startNodeName, virtualOutputStartType,weight);
                        std::tuple<std::string,std::string,double> undirectedEdgestartType(virtualInputStartType, startNodeName,weight);
                        std::tuple<std::string,std::string,double> edgeEndType(virtualInputEndType, endNodeName,weight);
                        std::tuple<std::string,std::string,double> undirectedEdgeEndType(endNodeName, virtualOutputEndType,weight);
                        // add the edge to the startType
                        if(!ret.first.contains(startType)){
                            ret.first[startType] = std::vector<std::tuple<std::string,std::string,double>>();
                        }
                        ret.first[startType].push_back(edgestartType);
                        if(undirectedTypeEdges){
                            ret.first[startType].push_back(undirectedEdgestartType);
                        }

                        // add the edge to the endType
                        if(!ret.first.contains(endType)){
                            ret.first[endType] = std::vector<std::tuple<std::string,std::string,double>>();
                        }
                        ret.first[endType].push_back(edgeEndType);
                        if(undirectedTypeEdges){
                            ret.first[endType].push_back(undirectedEdgeEndType);
                        }
                    } else {
                        //ignored because not in the subtypes
                    }
                    // add the edge with the contact times to the second vector in ret
                    ret.second.push_back(std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>(startNodeName, endNodeName, startType, endType, contactTimes, weight));
                    if(undirectedTypeEdges){
                        ret.second.push_back(std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>(endNodeName, startNodeName, endType, startType, contactTimes, weight));
                    }
                } else {
                    std::cout << "[ERROR] columns detected: " << entries.size() << " columns " <<std::endl;
                    throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: header doesn't have the right amount of columns(5 or 6 when considering interaction times) ");
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: file does not exists " + filename);
    }
    return ret;
}

std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>>> interactionContinuousContactsFileToEdgesListAndNodesByName(std::string filename, std::vector<std::string> subtypes,double maximumIntertypeTime, bool useEntrez, std::string granularity, std::unordered_map<std::string,std::vector<std::string>> typeToNodeNames , bool undirectedTypeEdges, double timestep){
    string line;
    std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>>> ret;
    // control if the granularity is valid
    if(granularity != "" && granularity != "type" && granularity != "node" && granularity != "typeAndNode"){
        throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: invalid granularity, it must be typeAndNode(finer) or type(coarser), or only node(no types)");
    }
    if(granularity == ""){
        granularity = "type";
    }

    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> entriesHeader = splitStringIntoVector(line, "\t");
            int indexTypeStart=-1, indexTypeEnd=-1, indexStartNode=-1, indexEndNode=-1, indexWeight=-1, indexContactTimes=-1;
            for(uint i = 0; i < entriesHeader.size(); i++){ //TODO change names in the header to be more general
                if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("starttype") != std::string::npos) {
                    indexTypeStart = i;
                }
                else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endtype") != std::string::npos) {
                    indexTypeEnd = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("startnodename") != std::string::npos) {
                    indexStartNode = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endnodename") != std::string::npos) {
                    indexEndNode = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("weight") != std::string::npos) {
                    indexWeight = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("contacttimes") != std::string::npos) {
                    indexContactTimes = i;
                }
                //startType	startNodeName	endType	endNodeName	weight
            }
            bool noContactTimes = false;
            if(indexTypeStart < 0 || indexTypeEnd < 0 || indexStartNode < 0 || indexEndNode < 0 || indexWeight < 0){
                throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: invalid file, the header does not contain a startType, or an endType, or a start node, or an end node, or a weight feature");
            }
            if(indexContactTimes < 0){
                noContactTimes = true;
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==5 || entries.size()==6){
                    std::string startNodeName,endNodeName;
                    if(!useEntrez){
                        startNodeName = entries[indexStartNode];
                        endNodeName = entries[indexEndNode];
                        
                    } else{
                        if(mapEnsembleToEntrez.contains(entries[indexStartNode]) && mapEnsembleToEntrez.contains(entries[indexEndNode])){
                            startNodeName = mapEnsembleToEntrez[entries[indexStartNode]];
                            endNodeName = mapEnsembleToEntrez[entries[indexEndNode]];
                        }
                    }
                    std::string startType = entries[indexTypeStart];
                    std::string endType = entries[indexTypeEnd];
                    if(typeToNodeNames.size() != 0 && vectorContains(subtypes, startType) && vectorContains(subtypes, endType)){
                        if(!vectorContains(typeToNodeNames[startType],startNodeName)){
                            std::cout << "[ERROR] start node <"<< startNodeName <<"> for type: " << startType << " is not in the specified network, aborting " <<std::endl;
                            throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: invalid file, the start node" + startNodeName + " is not in the type specified, aborting");
                        }

                        if(!vectorContains(typeToNodeNames[endType],endNodeName)){
                            std::cout << "[ERROR] end node <"<< endNodeName <<"> for type: " << endType << " is not in the specified network, aborting " <<std::endl;
                            throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: invalid file, the end node" + endNodeName + " is not in the type specified, aborting");
                        }
                    }

                    std::set<double> contactTimes;
                    double weight = std::stod( entries[indexWeight]);
                    if(noContactTimes){
                        // if no contact times are specified, then every time is a contact time, so all contact times from 0 to maximumIntertypeTime are added
                        for(int i = 0; i < maximumIntertypeTime; i++){
                            contactTimes.insert(i*timestep);
                        }
                    } else {
                        std::string contactTimesString = entries[indexContactTimes];
                        std::vector<std::string> splittedContactTimes = splitStringIntoVector(contactTimesString, ",");
                        for(auto iter = splittedContactTimes.cbegin(); iter!=splittedContactTimes.cend(); iter++){
                            double contactTime = std::stod(*iter);
                            // TODO change interTypeTime to be a double since it is a time, and I am now using the timestep to define the times and quantize the contacts
                            if(contactTime <= maximumIntertypeTime){
                                contactTimes.insert(contactTime);
                            } else {
                                std::cout << "[WARNING] contact time: " << contactTime << " is greater than the maximumIntertypeTime: " << maximumIntertypeTime << " ignoring it" <<std::endl;
                            }
                        }
                    }
                    // add the edge to the augmented graph, only if the startType and the endType are in the subtypes
                    if(vectorContains(subtypes, startType) && vectorContains(subtypes, endType)){
                        std::string virtualInputEndType = "";
                        std::string virtualOutputStartType = "";
                        // when undirectedTypeEdges is true, virtual input for the startType and virtual output for the endType are added to the respective graphs 
                        std::string virtualInputStartType = "";
                        std::string virtualOutputEndType = "";
                        if(granularity == "typeAndNode"){
                            virtualInputEndType = "v-in:" + startType + "_" + startNodeName;
                            virtualOutputStartType = "v-out:" + endType + "_" + endNodeName;
                            virtualInputStartType = "v-in:" + endType + "_" + endNodeName;
                            virtualOutputEndType = "v-out:" + startType + "_" + startNodeName;
                        } else if (granularity == "type"){
                            virtualInputEndType = "v-in:" + startType;
                            virtualOutputStartType = "v-out:" + endType;
                            virtualInputStartType = "v-in:" + endType;
                            virtualOutputEndType = "v-out:" + startType;
                        } else {
                            virtualInputEndType = "v-in:" + startNodeName;
                            virtualOutputStartType = "v-out:" + endNodeName;
                            virtualInputStartType = "v-in:" + endNodeName;
                            virtualOutputEndType = "v-out:" + startNodeName;
                        }
                        std::tuple<std::string,std::string,double> edgestartType(startNodeName, virtualOutputStartType,weight);
                        std::tuple<std::string,std::string,double> undirectedEdgestartType(virtualInputStartType, startNodeName,weight);
                        std::tuple<std::string,std::string,double> edgeEndType(virtualInputEndType, endNodeName,weight);
                        std::tuple<std::string,std::string,double> undirectedEdgeEndType(endNodeName, virtualOutputEndType,weight);
                        // add the edge to the startType
                        if(!ret.first.contains(startType)){
                            ret.first[startType] = std::vector<std::tuple<std::string,std::string,double>>();
                        }
                        ret.first[startType].push_back(edgestartType);
                        if(undirectedTypeEdges){
                            ret.first[startType].push_back(undirectedEdgestartType);
                        }

                        // add the edge to the endType
                        if(!ret.first.contains(endType)){
                            ret.first[endType] = std::vector<std::tuple<std::string,std::string,double>>();
                        }
                        ret.first[endType].push_back(edgeEndType);
                        if(undirectedTypeEdges){
                            ret.first[endType].push_back(undirectedEdgeEndType);
                        }
                    } else {
                        //ignored because not in the subtypes
                    }
                    // add the edge with the contact times to the second vector in ret
                    ret.second.push_back(std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>(startNodeName, endNodeName, startType, endType, contactTimes, weight));
                    if(undirectedTypeEdges){
                        ret.second.push_back(std::tuple<std::string, std::string, std::string, std::string, std::set<double>, double>(endNodeName, startNodeName, endType, startType, contactTimes, weight));
                    }
                } else {
                    std::cout << "[ERROR] columns detected: " << entries.size() << " columns " <<std::endl;
                    throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: header doesn't have the right amount of columns(5 or 6 when considering interaction times) ");
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::interactionContinuousContactsFileToEdgesListAndNodesByName: file does not exists " + filename);
    }

    return ret; 
}


std::map<std::string, std::string> getEnsembletoEntrezidMap(){
    string line;
    std::map<std::string,std::string> ret;
    std::string mapFilename = "resources/graphs/metapathwayNew/nodes.tsv";
    if(file_exists(mapFilename)){
        ifstream myfile (mapFilename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==4){
                    std::string Id = entries[0];
                    std::string Name = entries[1];
                    std::string Type = entries[2];
                    std::string Aliases = entries[3];
                    ret[Name] = Id;
                    
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::getEnsembletoEntrezidMap: file does not exists " + mapFilename);
    }
    return ret;
}

std::map<std::string, std::vector<std::string>> getFullNodesDescription(std::string filename){
    string line;
    // schema is #Id	Name	Type	Aliases
    std::map<std::string,std::vector<std::string>> ret;
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitStringIntoVector(line, "\t");
                if(entries.size()==4){
                    std::string Id = entries[0];
                    std::string Name = entries[1];
                    std::string Type = entries[2];
                    std::string Aliases = entries[3];
                    ret[Id] = std::vector<std::string>{Id,Name,Type,Aliases};
                    
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::getFullNodesDescription: file does not exists " + filename);
    }
    return ret;
}


std::vector<std::string> get_all(std::string const & root, std::string const & ext)
{
    std::vector<std::string> paths;
    for (auto &p : std::filesystem::recursive_directory_iterator(root))
    {
        if (p.path().extension() == ext)
            paths.push_back(root + "/" + p.path().stem().string() + ext);
    }
    return paths;
} 

void saveNodeValues(std::string folderName, int iteration, std::string cellName, std::vector<double> nodeValues,std::vector<std::string> nodeNames, bool useEntrez, std::string nodesDescriptionFile){
    std::string outputFilename = folderName + "/" + cellName + "--"+std::to_string(iteration) + ".tsv";
    std::ofstream outfile(outputFilename,ios::out|ios::trunc);

    if (!outfile.is_open()) {
        std::cout << "Unable to open file " << outputFilename << std::endl;
        return;
    }
    if(nodesDescriptionFile.length()==0){
        if(useEntrez)
            nodesDescriptionFile = "resources/graphs/metapathwayNew/nodes.tsv";
        else  // no nodes description file is used
            nodesDescriptionFile = "";
    } else {
        if(!file_exists(nodesDescriptionFile)){
            throw std::invalid_argument("utilities::saveNodeValues: file does not exists " + nodesDescriptionFile);
        }
    }
    //auto mapToEnsemble = getEnsembletoEntrezidMap();
    //header
    outfile << "nodeID\tnodeName\ttype\talias\tnodeValue\n";
    //body
    if(useEntrez || nodesDescriptionFile.length()!=0){
        auto mapToEverything = getFullNodesDescription(nodesDescriptionFile);
        for(uint i = 0; i < nodeValues.size(); i++){
            if(mapToEverything.contains(nodeNames[i]))
                outfile<<mapToEverything.at(nodeNames[i])[0]<<"\t"<<mapToEverything.at(nodeNames[i])[1]<<"\t"<<mapToEverything.at(nodeNames[i])[2]<<"\t"<<mapToEverything.at(nodeNames[i])[3]<<"\t"<< std::to_string(nodeValues[i]);
            else {
                auto splittedVirtual = splitStringIntoVector(nodeNames[i], ":");
                if(splittedVirtual[0]=="v-in"){
                    outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"virtual-input\t"<<splittedVirtual[1]<<"\t"<<std::to_string(nodeValues[i]);
                } else if(splittedVirtual[0]=="v-out"){
                    outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"virtual-output\t"<<splittedVirtual[1]<<"\t"<<std::to_string(nodeValues[i]);
                } else{ //when the node names are not genes but something else
                    outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"nodes in the graph\t"<<nodeNames[i]<<"\t"<<std::to_string(nodeValues[i]);
                }
            }
            outfile << std::endl;
        }
    }  
    else{
        for(uint i = 0; i < nodeValues.size(); i++){
            auto splittedVirtual = splitStringIntoVector(nodeNames[i], ":");
            if(splittedVirtual[0]=="v-in"){
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"virtual-input\t"<<splittedVirtual[1]<<"\t"<<std::to_string(nodeValues[i]);
            } else if(splittedVirtual[0]=="v-out"){
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"virtual-output\t"<<splittedVirtual[1]<<"\t"<<std::to_string(nodeValues[i]);
            } else{ //when the node names are not genes but something else
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"nodes in the graph\t"<<nodeNames[i]<<"\t"<<std::to_string(nodeValues[i]);
            }
            outfile << std::endl;
        }
    }

    // for (const auto& row : data) {
    //     for (const auto& item : row) {
    //         outfile << item << ",";
    //     }
    //     outfile << "\n";
    // }

    outfile.close();
}

void saveNodeValues(std::string folderName, int iterationOuter, int intraIteration, std::string cellName, std::vector<double> nodeValues,std::vector<std::string> nodeNames, bool useEntrez, std::string nodesDescriptionFile){
    std::string outputFilename = folderName + "/" + cellName + "--"+std::to_string(iterationOuter + intraIteration) + ".tsv";
    std::ofstream outfile(outputFilename,ios::out|ios::trunc);

    if (!outfile.is_open()) {
        std::cout << "Unable to open file " << outputFilename << std::endl;
        throw std::invalid_argument("utilities::saveNodeValues: unable to open output file " + outputFilename);
    }

    if(nodesDescriptionFile.length()==0 && useEntrez){
            nodesDescriptionFile = "resources/graphs/metapathwayNew/nodes.tsv";
    } else if(nodesDescriptionFile.length()!=0 && !file_exists(nodesDescriptionFile)){
        throw std::invalid_argument("utilities::saveNodeValues: file does not exists " + nodesDescriptionFile);
    }

    std::map<std::string, std::vector<std::string>> mapToEverything;
    if(useEntrez || nodesDescriptionFile.length()!=0){
        mapToEverything = getFullNodesDescription(nodesDescriptionFile);
    } else {
        mapToEverything = std::map<std::string, std::vector<std::string>>();
    }

    //header
    outfile << "nodeID\tnodeName\ttype\talias\tnodeValue\n";
    //body
    for(uint i = 0; i < nodeValues.size(); i++){
        if(mapToEverything.size() !=0 && mapToEverything.contains(nodeNames[i]))
            outfile<<mapToEverything.at(nodeNames[i])[0]<<"\t"<<mapToEverything.at(nodeNames[i])[1]<<"\t"<<mapToEverything.at(nodeNames[i])[2]<<"\t"<<mapToEverything.at(nodeNames[i])[3]<<"\t"<< std::to_string(nodeValues[i]);
        else {
            auto splittedVirtual = splitStringIntoVector(nodeNames[i], ":");
            if(splittedVirtual[0]=="v-in"){
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"virtual-input\t"<<splittedVirtual[1]<<"\t"<<std::to_string(nodeValues[i]);
            } else if(splittedVirtual[0]=="v-out"){
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"virtual-output\t"<<splittedVirtual[1]<<"\t"<<std::to_string(nodeValues[i]);
            } else{ //when the node names are not genes but something else
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"nodes in the graph\t"<<nodeNames[i]<<"\t"<<std::to_string(nodeValues[i]);
            }
        }
        outfile << std::endl;
    }
}

void saveNodeValuesWithTime(std::string folderName,int iterationOuter, int intraIteration, std::string cellName, std::vector<double> nodeValues,std::vector<std::string> nodeNames, bool useEntrez, std::string nodesDescriptionFile, double timestep){
    std::string outputFilename = folderName + "/" + cellName + "--"+std::to_string(iterationOuter + intraIteration) + ".tsv";
    std::ofstream outfile(outputFilename,ios::out|ios::trunc);

    if (!outfile.is_open()) {
        std::cout << "Unable to open file " << outputFilename << std::endl;
        throw std::invalid_argument("utilities::saveNodeValues: unable to open output file " + outputFilename);
    }

    //TODO specifiy if the nodes description are in a folder with the graphs description files?
    if(nodesDescriptionFile.length()==0 && useEntrez){
            nodesDescriptionFile = "resources/graphs/metapathwayNew/nodes.tsv";
    } else if(nodesDescriptionFile.length()!=0 && !file_exists(nodesDescriptionFile)){
        throw std::invalid_argument("utilities::saveNodeValues: file does not exists " + nodesDescriptionFile);
    }

    std::map<std::string, std::vector<std::string>> mapToEverything;
    if(useEntrez || nodesDescriptionFile.length()!=0){
        mapToEverything = getFullNodesDescription(nodesDescriptionFile);
    } else {
        mapToEverything = std::map<std::string, std::vector<std::string>>();
    }

    //header
    outfile << "nodeID\tnodeName\ttype\talias\tnodeValue\ttime\n";
    //body
    for(uint i = 0; i < nodeValues.size(); i++){
        if(mapToEverything.size() !=0 && mapToEverything.contains(nodeNames[i]))
            outfile<<mapToEverything.at(nodeNames[i])[0]<<"\t"<<mapToEverything.at(nodeNames[i])[1]<<"\t"<<mapToEverything.at(nodeNames[i])[2]<<"\t"<<mapToEverything.at(nodeNames[i])[3]<<"\t"<< std::to_string(nodeValues[i]) << "\t" << std::to_string((iterationOuter + intraIteration) *timestep);
        else {
            auto splittedVirtual = splitStringIntoVector(nodeNames[i], ":");
            if(splittedVirtual[0]=="v-in"){
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"virtual-input\t"<<splittedVirtual[1]<<"\t"<<std::to_string(nodeValues[i])<< "\t" << std::to_string((iterationOuter + intraIteration) *timestep);
            } else if(splittedVirtual[0]=="v-out"){
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"virtual-output\t"<<splittedVirtual[1]<<"\t"<<std::to_string(nodeValues[i])<< "\t" << std::to_string((iterationOuter + intraIteration) *timestep);
            } else{ //when the node names are not genes but something else
                outfile<<nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<"nodes in the graph\t"<<nodeNames[i]<<"\t"<<std::to_string(nodeValues[i])<< "\t" << std::to_string((iterationOuter + intraIteration) *timestep);
            }
        }
        outfile << std::endl;
    }
}


double vectorNorm(std::vector<double> vec){
    double norm=0;
    for (uint i = 0; i < vec.size(); ++i) {
        norm+=vec[i]*vec[i];
    }
    norm=sqrt(norm);
    return norm;
}