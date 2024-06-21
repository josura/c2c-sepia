#pragma once
#include <string>
#include <vector>
#include "Computation.h"

class Checkpoint {
public:
    Checkpoint();
    ~Checkpoint();

    /**
     * \brief save the state of the computation of a type in a file
     * \param type: the type of the computation
     * \param interIteration: the number of the interation
     * \param intraIteration: the number of the intra iteration
     * \param currentComputation: the computation to save
     */
    void saveState(const std::string type, const int interIteration, const int intraIteration, const Computation* currentComputation);
    /**
     * \brief remove all the checkpoints of a given type, deleting the checkpoint files
     */
    void cleanCheckpoints(const std::string type);
    /**
     * \brief load the state of the computation of a type from a file
     * \param type: the type of the computation
     * \param interIteration: the number of the interation
     * \param intraIteration: the number of the intra iteration
     * \param computation: the computation to load
     */
    void loadState(const std::string type, int& interIteration, int& intraIteration, Computation* computation);

    // Add your variables and methods here
    void setCheckPointFolder(const std::string& folder);

private:
    std::string checkPointFolder;
};
