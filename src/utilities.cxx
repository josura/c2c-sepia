#include "utilities.h"
#include <algorithm>
#include <boost/token_functions.hpp>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iterator>
#include <map>
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
    val =  randomNumber(-INTMAX, INTMAX);
}
void setRandom(double& val) { 
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


std::vector<std::string> splitString(std::string toSplit , std::string delimiter){

    vector<string> tokens;
    boost::split(tokens, toSplit, boost::is_any_of(delimiter));
    return tokens;
}


boost::tokenizer<boost::char_delimiters_separator<char>> splitStringTokenizer(std::string toSplit , char delimiter){
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
                std::vector<std::string> entries = splitString(line, "\t");
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
            std::vector<std::string> entriesHeader = splitString(line, "\t");
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
                throw std::invalid_argument("invalid file, the header does not contain a start, an end and a weight feature");
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitString(line, "\t");
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
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::edgeFileEdgesListByIndex: file does not exists " + filename);
    }
    return std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>> (nameRet,ret);
}


std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> logFoldChangeMatrixToCellVectors(std::string filename, std::vector<std::string> finalNames,bool useEntrez){
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
            std::vector<std::string> splittedHeader = splitString(line, "\t");  //could already be used as the cellnames vector
            for (int i = 0; i < SizeToInt( splittedHeader.size()); i++) {
                cellNames.push_back(splittedHeader[i]);
                ret.push_back(std::vector<double>(finalNames.size(),0));
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitString(line, "\t");
                if(entries.size()>1){
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
                }
            }
            myfile.close();
            std::cout << "[LOG] No node in the metapathway for genes: " << std::endl;
            for(auto iter = discardedGenes.cbegin();iter!=discardedGenes.cend();iter++){
                std::cout << "," << *iter;
            }
            std::cout << std::endl <<"[LOG] discarding logfold for the genes not in the metapathway" << std::endl;
            
        }
    } else {
        throw std::invalid_argument("utilities::edgeFileEdgesListByIndex: file does not exists " + filename);
    }
    return std::tuple<std::vector<std::string>,std::vector<std::string>,std::vector<std::vector<double>>> (geneNames,cellNames,ret);

}


//TODO, understand what file or files(maybe a directory) should be read into the program, dependent on how the cells are represented
//TODO, understand if the translation from ensemble gene names to entrez should be done here
//TODO, filtering genes also since I have seen nodes not in the metapathway
std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> cellInteractionFileToEdgesListAndNodesByName(std::string filename,bool useEntrez){
    string line;
    std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>> ret;
    auto mapEnsembleToEntrez = getEnsembletoEntrezidMap();
    if(file_exists(filename)){
        ifstream myfile (filename);
        if (myfile.is_open())
        {
            getline (myfile,line);  // first line is header IMPORTANT
            std::vector<std::string> entriesHeader = splitString(line, "\t");
            int indexCellStart=-1,indexCellEnd=-1,indexLigandStart=-1,indexReceptorEnd=-1,indexWeight=-1;
            for(uint i = 0; i < entriesHeader.size(); i++){
                if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("startcell") != std::string::npos) {
                    indexCellStart = i;
                }
                else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("endcell") != std::string::npos) {
                    indexCellEnd = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("geneligand") != std::string::npos) {
                    indexLigandStart = i;
                }else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("genereceptor") != std::string::npos) {
                    indexReceptorEnd = i;
                } else if (boost::algorithm::to_lower_copy(entriesHeader[i]).find("weight") != std::string::npos) {
                    indexWeight = i;
                }
            }
            if(indexCellStart < 0 || indexCellEnd < 0 || indexLigandStart < 0 || indexReceptorEnd < 0 || indexWeight < 0){
                throw std::invalid_argument("invalid file, the header does not contain a startcell, or an endcell, or a Ligand gene, or a receptor gene, or a weight feature");
            }
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitString(line, "\t");
                if(entries.size()==5){
                    std::string geneLigand,geneReceptor;
                    if(!useEntrez){
                        geneLigand = entries[indexLigandStart];
                        geneReceptor = entries[indexReceptorEnd];    

                        std::string startCell = entries[indexCellStart];
                        std::string endCell = entries[indexCellEnd];
                        double weight = std::stod( entries[indexWeight]);
                        std::string virtualInputEndCell = "v-in:" + endCell;
                        std::string virtualOutputStartCell = "v-out:" + startCell;
                        std::tuple<std::string,std::string,double> edgeStartCell(geneLigand, virtualOutputStartCell,weight);
                        std::tuple<std::string,std::string,double> edgeEndCell(virtualInputEndCell, geneReceptor,weight);
                        if(ret.contains(startCell)){
                            ret[startCell].push_back(edgeStartCell);
                        }else{
                            ret[startCell] = std::vector<std::tuple<std::string,std::string,double>>();
                            ret[startCell].push_back(edgeStartCell);
                        }

                        if(ret.contains(endCell)){
                            ret[endCell].push_back(edgeEndCell);
                        }else{
                            ret[endCell] = std::vector<std::tuple<std::string,std::string,double>>();
                            ret[endCell].push_back(edgeEndCell);
                        }
                    } else{
                        if(mapEnsembleToEntrez.contains(entries[indexLigandStart]) && mapEnsembleToEntrez.contains(entries[indexReceptorEnd])){
                            geneLigand = mapEnsembleToEntrez[entries[indexLigandStart]];
                            geneReceptor = mapEnsembleToEntrez[entries[indexReceptorEnd]];

                            std::string startCell = entries[indexCellStart];
                            std::string endCell = entries[indexCellEnd];
                            double weight = std::stod( entries[indexWeight]);
                            std::string virtualInputEndCell = "v-in:" + endCell;
                            std::string virtualOutputStartCell = "v-out:" + startCell;
                            std::tuple<std::string,std::string,double> edgeStartCell(geneLigand, virtualOutputStartCell,weight);
                            std::tuple<std::string,std::string,double> edgeEndCell(virtualInputEndCell, geneReceptor,weight);
                            if(ret.contains(startCell)){
                                ret[startCell].push_back(edgeStartCell);
                            }else{
                                ret[startCell] = std::vector<std::tuple<std::string,std::string,double>>();
                                ret[startCell].push_back(edgeStartCell);
                            }

                            if(ret.contains(endCell)){
                                ret[endCell].push_back(edgeEndCell);
                            }else{
                                ret[endCell] = std::vector<std::tuple<std::string,std::string,double>>();
                                ret[endCell].push_back(edgeEndCell);
                            }
                        }
                    }
                }
            }
            myfile.close();
        }
    } else {
        throw std::invalid_argument("utilities::cellInteractionFileToEdgesListAndNodesByName: file does not exists " + filename);
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
                std::vector<std::string> entries = splitString(line, "\t");
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