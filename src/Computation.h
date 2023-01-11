

#include "Matrix.h"
#include <map>
#include <string>
#include <vector>

class Computation{
    private:
        std::vector<double> input,output;
        Matrix<double> W;
        Matrix<double> Wstar;
        std::map<std::string, double> geneMapToNode;
        std::vector<std::string> cellTypes;
    public:
        Computation();
        Computation(std::vector<double> _input);
        Computation(std::vector<double> _input, Matrix<double> _W);
        Computation(std::vector<double> _input, Matrix<double> _W,std::vector<std::string> _cellTypes);
};