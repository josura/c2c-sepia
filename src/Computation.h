

#include "Matrix.h"
#include <map>
#include <string>
#include <tuple>
#include <vector>

class Computation{
    private:
        std::vector<double> input,output;
        Matrix<double> W;
        Matrix<double> Wstar;
        std::map<std::string, int> geneMapToNode;
        std::vector<std::string> cellTypes;
    public:
        Computation();
        Computation(std::vector<double>& _input);
        Computation(std::vector<double>& _input, Matrix<double>& _W);
        Computation(std::vector<double>& _input, Matrix<double>& _W,std::vector<std::string>& _cellTypes);
        void augmentW(std::vector<std::tuple<std::string,std::string,double>>);
        std::vector<double> computePerturbation();
        std::vector<double> computeAugmentedPerturbation(); //taking into account virtual nodes in the augmented metapathway

        
};