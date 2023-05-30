#pragma once

#include "DissipationModel.h"
#include "Matrix.h"
#include "WeightedEdgeGraph.h"
#include <map>
#include <string>
#include <tuple>
#include <vector>

class Computation{
    private:
        std::vector<double> input,output,inputAugmented,outputAugmented;
        WeightedEdgeGraph* metapathway;
        WeightedEdgeGraph* augmentedMetapathway;
        std::vector<std::string> cellTypes;
        std::string localCellType;
        bool armaInitializedNotAugmented = false, armaInitializedAugmented = false;
        arma::Col<double> InputArma;
        arma::Mat<double> pseudoInverseArma;
        arma::Col<double> InputAugmentedArma;
        arma::Mat<double> pseudoInverseAugmentedArma;
        std::map<std::string,int> nodeToIndex;
        DissipationModel* dissipationModel=nullptr;
    public:
        Computation();
        ~Computation();
        Computation(std::string _thisCellType, const std::vector<double>& _input);   // useless???

        /*
        Computation without knowledge of the other cell types, this part can be seen as the classical algorithm without additional computation for message passing
        between cells, only intra-cell propagation
        @param std::string _thisCellType: the celltype of this computation, this information will be used as the unique name for the Agent
        @param const std::vector<double>& _input: input vector of the nodes values, initially the one passed in the input
        @param const Matrix<double>& _W: the adjacency matrix along the values of every edge in the graph that it represents
        @param const std::vector<std::string>& metapathwayNames: the graph nodes names, in order defined by the adjacency matrix
        */
        Computation(std::string _thisCellType,const std::vector<double>& _input, const Matrix<double>& _W, const std::vector<std::string>& metapathwayNames);


        /*
        Computation without knowledge of the other cell types, this part can be seen as the classical algorithm without additional computation for message passing
        between cells, only intra-cell propagation
        @param std::string _thisCellType: the celltype of this computation, this information will be used as the unique name for the Agent
        @param const std::vector<double>& _input: input vector of the nodes values, initially the one passed in the input
        @param const Matrix<double>& _W: the adjacency matrix along the values of every edge in the graph that it represents
        @param const std::vector<std::string>& metapathwayNames: the graph nodes names, in order defined by the adjacency matrix
        */
        Computation(std::string _thisCellType,const std::vector<double>& _input, WeightedEdgeGraph* _metapathway, const std::vector<std::string>& metapathwayNames);
        
        /*
        Augment the metapathway with celltypes and a new set of edges from virtual nodes in the augmented metapathway to the metapathway(virtual inputs and virtual outputs) 
        @param const std::vector<std::string>& _cellTypes: the celltypes other than this celltype, the other agents in the network
        @param 
        */
        void augmentMetapathway(const std::vector<std::string>&,const std::vector<std::pair<std::string,std::string>>& newEdgesList =std::vector<std::pair<std::string,std::string>>(), const std::vector<double>& newEdgesValue = std::vector<double>(), bool includeSelfVirtual=false);
        void augmentMetapathwayNoComputeInverse(const std::vector<std::string>&,const std::vector<std::pair<std::string,std::string>>& newEdgesList =std::vector<std::pair<std::string,std::string>>(), const std::vector<double>& newEdgesValue = std::vector<double>(), bool includeSelfVirtual=false);
        void addEdges(const std::vector<std::pair<std::string,std::string>>& , const std::vector<double>& , bool bothDirections = false);
        void addEdges(const std::vector<std::tuple<std::string,std::string,double>>&  , bool bothDirections = false);
        void addEdges(const std::vector<std::pair<int,int>>& , const std::vector<double>& , bool bothDirections = false);
        void addEdges(const std::vector<std::tuple<int,int,double>>&  , bool bothDirections = false);
        std::vector<double> computePerturbation();
        std::vector<double> computeAugmentedPerturbation(); //taking into account virtual nodes in the augmented metapathway
        std::vector<double> computeAugmentedPerturbationNorm(); //taking into account normalization
        std::vector<double> computeAugmentedPerturbationDissipatedPow2(); //taking into account dissipation after every iteration(a kind of starvation probably) subtracting the pow2 of every perturbation before the computation
        std::vector<double> computeAugmentedPerturbationDissipatedAfterCompute(double timeStep); //taking into account dissipation after every iteration(Dissipation model), dissipation after the computation of the perturbated value
        std::vector<double> computeAugmentedPerturbationDissipatedBeforeCompute(double timeStep); //taking into account dissipation after every iteration (Dissipation model), dissipation before the computation of the perturbated value
        std::vector<double> computeAugmentedPerturbationSaturated(const std::vector<double>& saturationsVector = std::vector<double>()); //taking into account saturation(hyperbolic tangent and scaling) after every iteration
        std::vector<double> computeAugmentedPerturbationSaturatedAndDissipatedBeforeCompute(double timeStep,const std::vector<double>& saturationsVector = std::vector<double>()); //taking into account saturation(hyperbolic tangent and scaling) and dissipation after every iteration
        std::vector<double> computeAugmentedPerturbationEnhanced1(const std::vector<double>& saturationsVector = std::vector<double>()); //taking into account saturation(hyperbolic tangent and scaling) and dissipation after every iteration
        std::vector<double> computeAugmentedPerturbationEnhanced2(const std::vector<double>& saturationsVector = std::vector<double>()); //taking into account saturation(hyperbolic tangent and scaling) and dissipation after every iteration
        std::pair<std::string,double> getMapVirtualOutputsToCellInputs(); //TODO
        void updateInput(const std::vector<double>& newInp = std::vector<double>(), bool augmented = false);

        // get sets
        
        std::vector<double> getInput()const{return input;}
        std::vector<double> getOutput()const{return output;}
        std::vector<double> getInputAugmented()const{return inputAugmented;}
        std::vector<double> getOutputAugmented()const{return outputAugmented;}
        WeightedEdgeGraph* getMetapathway()const{return metapathway;}
        WeightedEdgeGraph* getAugmentedMetapathway()const{return augmentedMetapathway;}
        std::vector<std::string> getCellTypes()const{return cellTypes;}
        std::string getLocalCellType()const{return localCellType;}
        bool isInitializedArmaNotAugmented()const{return armaInitializedNotAugmented;}
        bool isInitializedArmaAugmented()const{return armaInitializedAugmented;}
        arma::Col<double> getInputArma()const{return InputArma;}
        arma::Mat<double> getPseudoInverseArma()const{return pseudoInverseArma;}
        arma::Col<double> getInputAugmentedArma()const{return InputAugmentedArma;}
        arma::Mat<double> getPseudoInverseAugmentedArma()const{return pseudoInverseAugmentedArma;}

        double getVirtualInputForCell(std::string celltype)const;
        double getVirtualOutputForCell(std::string celltype)const;
        void setInputVinForCell(std::string celltype, double value);
        void setInputVoutForCell(std::string celltype, double value);
        void setDissipationModel(DissipationModel* dissipationModel);

        //optimization
        void freeAugmentedGraphs();

        // operators
        Computation& operator=( const Computation& );
        Computation copy()const;
        void assign(const Computation&);
        
};