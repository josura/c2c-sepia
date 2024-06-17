#ifndef CHECKPOINT_HXX
#define CHECKPOINT_HXX

#include <string>

class Checkpoint {
public:
    Checkpoint();
    ~Checkpoint();

    void saveState(const std::string& filename);
    void loadState(const std::string& filename);

    // Add your variables and methods here

private:
    // Add your private variables and methods here
};

#endif // CHECKPOINT_HXX