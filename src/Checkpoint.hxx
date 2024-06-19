#pragma once
#include <string>
#include <vector>
#include "Computation.h"

class Checkpoint {
public:
    Checkpoint();
    ~Checkpoint();

    void saveState(const std::string type, const int interIteration, const int intraIteration, const Computation* currentComputation);
    void loadState(const std::string type, int& interIteration, int& intraIteration, Computation* computation);

    // Add your variables and methods here
    void setCheckPointFolder(const std::string& folder);

private:
    std::string checkPointFolder;
};
