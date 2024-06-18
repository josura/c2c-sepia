#pragma once
#include <string>

class Checkpoint {
public:
    Checkpoint();
    ~Checkpoint();

    void saveState(const std::string& filename);
    void loadState(const std::string& filename);

    // Add your variables and methods here

private:
    std::string checkPointFolder;
};
