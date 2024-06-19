#include<Checkpoint.hxx>
#include<Computation.h>
#include<string>
#include<utilities.h>

Checkpoint::Checkpoint() {
    this->checkPointFolder = "checkpoints/";
    // control if the folder exists
    // if not, create it
    if (folderExists(this->checkPointFolder) == false)
    {
        createFolder(this->checkPointFolder);
    }
    
}

Checkpoint::~Checkpoint() {
    // Add your code here
}

void Checkpoint::saveState(const std::string type, const int interIteration, const int intraIteration, const Computation* currentComputation) {
    std::string fileName = this->checkPointFolder + "checkpoint_" + type + "_" + std::to_string(interIteration) + "_" + std::to_string(intraIteration) + ".tsv";
    std::ofstream file(fileName);
    std::vector<std::string> nodeNames = currentComputation->getAugmentedGraph()->getNodeNames();
    std::vector<double> nodeValues = currentComputation->getOutputAugmented();
    if (file.is_open())
    {
        //header
        file << "interIteration\tintraIteration\tnodeName\tnodeValue\n";
        //body
        for(uint i = 0; i < nodeValues.size(); i++){
            file<<interIteration<< "\t"<< intraIteration<<"\t"<< nodeNames[i]<<"\t"<<nodeNames[i]<<"\t"<<std::to_string(nodeValues[i]);
            file << std::endl;
        }
        file.close();
    }
    else
    {
        std::cerr << "[ERROR] Checkpoint::saveState: Unable to open file " << fileName << std::endl;
    }
}

void Checkpoint::loadState(const std::string type, int& interIteration, int& intraIteration, Computation* computation) {
    std::string fileName = this->checkPointFolder + "checkpoint_" + type + "_" + std::to_string(interIteration) + "_" + std::to_string(intraIteration) + ".tsv";
    std::ifstream file(fileName);
    if (file.is_open())
    {
        std::string line;
        std::getline(file, line); // skip header
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::string interIterationStr;
            std::string intraIterationStr;
            std::string nodeName;
            std::string nodeValueStr;
            iss >> interIterationStr >> intraIterationStr >> nodeName >> nodeValueStr;
            double nodeValue = std::stod(nodeValueStr);
            interIteration = std::stoi(interIterationStr);
            intraIteration = std::stoi(intraIterationStr);
            computation->setInputNodeValue(nodeName, nodeValue);
        }
        file.close();
    }
    else
    {
        std::cerr << "[ERROR] Checkpoint::loadState: Unable to open file " << fileName << std::endl;
    }
}