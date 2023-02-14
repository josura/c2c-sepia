#include "utilities.h"
#include <algorithm>
#include <boost/token_functions.hpp>
#include <cstddef>
#include <cstdlib>
#include <iostream>
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

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
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
            while ( getline (myfile,line) )
            {
                std::vector<std::string> entries = splitString(line, "\t");
                if(entries.size()==3){
                    std::string node1 = entries[0];
                    std::string node2 = entries[1];
                    double weight = std::stod( entries[2]);
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
