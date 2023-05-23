#include <boost/program_options/value_semantic.hpp>
#include <iostream>
#include <boost/program_options.hpp>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include "Computation.h"
#include "DissipationModel.h"
#include "DissipationModelPow.h"
#include "DissipationModelRandom.h"
#include "DissipationModelScaled.h"
//#include "DissipationModelPeriodic.h"
#include "WeightedEdgeGraph.h"
#include "utilities.h"


void printHelp(){
    //TODO fix this help
    std::cout << "usage: ./c2c-sepia --fMETAPATHWAY <metapathway>.tsv --fLogfoldPerCell <logfoldPerCell>.tsv [<subcelltypes>.txt] --dirCellInteraction <cellInteractionFolder>(containing .tsv files)]\nFILE STRUCTURE SCHEMA:\nmetapathway.tsv\nstart\tend\tweight\n<gene1>\t<gene2>\t <0.something>\n...\n\n\nlogfoldPerCell.tsv\n\tcell1\tcell2\t...\tcellN\ngene1\t<lfc_cell1:gene1>\t<lfc_cell2:gene1>\t...\t<lfc_cellN:gene1>\ngene1\t<lfc_cell1:gene2>\t<lfc_cell2:gene2>\t...\t<lfc_cellN:gene2>\n...\n\n\ncelltypesInteraction.tsv\nstartCell:geneLigand\tendCell:geneReceptor\tweight\n<cell1:geneLigand>\t<cell2:genereceptor>\t <0.something>\n...\n\n\nsubcelltypes.txt\ncell1\ncell3\n..."<<std::endl;
    std::cout << "LEGEND:\n <> := placeholder for the name of the file\n[] := optional\n{} := at least one"<<std::endl;
}

int main(int argc, char** argv ) {
    //program options
    bool ensembleGeneNames=false;
    bool sameCellCommunication=false;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "print help section")//<logfoldPerCell>.tsv [<subcelltypes>.txt] [<celltypesInteraction>.tsv]\nFILE STRUCTURE SCHEMA:\nmetapathway.tsv\nstart end weight\n<gene1> <gene2>  <0.something>\n...\n\n\nlogfoldPerCell.tsv\n cell1 cell2 ... cellN\ngene1 <lfc_cell1:gene1> <lfc_cell2:gene1> ... <lfc_cellN:gene1>\ngene1 <lfc_cell1:gene2> <lfc_cell2:gene2> ... <lfc_cellN:gene2>\n...\n\n\ncelltypesInteraction.tsv\nstartCell:geneLigand endCell:geneReceptor weight\n<cell1:geneLigand> <cell2:genereceptor>  <0.something>\n...\n\n\nsubcelltypes.txt\ncell1\ncell3\n...")
        ("fMETAPATHWAY", po::value<std::string>()->required(), "metapathway filename, for an example metapathway see in resources")
        ("fLogfoldPerCell", po::value<std::string>()->required(), "logfoldPerCell matrix filename, for an example see in data")
        ("dirCellInteraction", po::value<std::string>(), "directory for the cell interactions, for an example see in data")
        ("ensembleGeneNames",po::bool_switch(&ensembleGeneNames),"use ensemble gene names, since the metapathway used in resources have entrez_ids, a map will be done from ensemble to entrez, the map is available in resources")
        ("sameCellCommunication",po::bool_switch(&sameCellCommunication),"use same cell communication, since it is not permitted as the standard definition of the model, this adds a virtual node for the same cell type")
        ("output",po::value<std::string>()->required(),"output folder for output of the algorithm at each iteration")
        ("intercellIterations",po::value<uint>(),"number of iterations for intercell communication")
        ("intracellIterations",po::value<uint>(),"number of iterations for intracell communication")
        ("dissipationModel",po::value<std::string>(),"the dissipation model for the computation, available models are: 'none (default)','power','random','periodic','scaled'")
        ("dissipationModelParameters",po::value<std::vector<double>>()->multitoken(),"the parameters for the dissipation model, for the power dissipation indicate the base, for the random dissipation indicate the min and max value, for the periodic dissipation indicate the period")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    std::string filename,celltypesFilename,celltypesInteractionFoldername,cellLogFoldMatrixFilename,outputFoldername;
    uint intercellIterations,intracellIterations;
    DissipationModel* dissipationModel = nullptr;

    if (vm.count("help")) {
        //printHelp();
        std::cout << desc << std::endl;
        return 1;
    }

    if (vm.count("intercellIterations")) {
        std::cout << "[LOG] iterations intercell set to " 
    << vm["intercellIterations"].as<std::string>() << ".\n";
        intercellIterations = vm["intercellIterations"].as<uint>();
    } else {
        std::cout << "[LOG] iterations intercell not set, set to default: 10 iterations \n";
        intercellIterations = 10;
    }
    if (vm.count("intracellIterations")) {
        std::cout << "[LOG] iterations intracell set to " 
    << vm["intracellIterations"].as<std::string>() << ".\n";
        intracellIterations = vm["intracellIterations"].as<uint>();
    } else {
        std::cout << "[LOG] iterations intracell not set, set to default: 5 iterations \n";
        intracellIterations = 5;
    }

    if (vm.count("fMETAPATHWAY")) {
        std::cout << "[LOG] file for the metapathway was set to " 
    << vm["fMETAPATHWAY"].as<std::string>() << ".\n";
        filename = vm["fMETAPATHWAY"].as<std::string>();
        if(!fileExistsPath(filename)){
            std::cerr << "[ERROR] file for the metapathway do not exist: aborting"<<std::endl;
            return 1;
        }
    } else {
        std::cout << "[ERROR] fMETAPATHWAY file was not set. Aborting\n";
        return 1;
    }
    if (vm.count("fLogfoldPerCell")) {
        std::cout << "[LOG] file for the logfoldPerCell matrix was set to " 
    << vm["fLogfoldPerCell"].as<std::string>() << ".\n";
        cellLogFoldMatrixFilename = vm["fLogfoldPerCell"].as<std::string>();
        if(!fileExistsPath(cellLogFoldMatrixFilename)){
            std::cerr << "[ERROR] file for the logfoldPerCell does not exist: aborting"<<std::endl;
            return 1;
        }
    } else {
        std::cerr << "[ERROR] fLogfoldPerCell file was not set. Aborting\n";
        return 1;
    }
    if (vm.count("dirCellInteraction")) {
        std::cout << "[LOG] folder for the cell interactions was set to " 
    << vm["dirCellInteraction"].as<std::string>() << ".\n";
        celltypesInteractionFoldername = vm["dirCellInteraction"].as<std::string>();
        if(!folderExists(celltypesInteractionFoldername)){
            std::cerr << "[ERROR] folder for the cell interactions do not exist: aborting"<<std::endl;
            return 1;
        }
    } else {
        std::cout << "[LOG] dirCellInteraction folder was not set. computing without taking into account cell interactions\n";
        //TODO
    }
    if (vm.count("output")) {
        std::cout << "[LOG] output folder  was set to " 
    << vm["output"].as<std::string>() << ".\n";
        outputFoldername = vm["output"].as<std::string>();
        if(!folderExists(outputFoldername)){
            std::cerr << "[ERROR] folder for the output do not exist: aborting"<<std::endl;
            return 1;
        }
    } else {
        std::cout << "[LOG] output folder was not set. aborting\n";
        return 1;
        //TODO
    }
    if (vm.count("dissipationModel")) {
        std::cout << "[LOG] dissipation model was set to "
    << vm["dissipationModel"].as<std::string>() << ".\n";
        std::string dissipationModelName = vm["dissipationModel"].as<std::string>();
        if(dissipationModelName == "none"){
            std::cout << "[LOG] dissipation model set to default (none)\n";
            dissipationModel = new DissipationModel();
        } else if(dissipationModelName == "power"){
            if (vm.count("dissipationModelParameters")) {
                std::cout << "[LOG] dissipation model parameters for power dissipation were declared to be" << vm["dissipationModelParameters"].as<std::vector<double>>()[0] << ".\n";
                std::vector<double> dissipationModelParameters = vm["dissipationModelParameters"].as<std::vector<double>>();
                if(dissipationModelParameters.size() == 1){
                    dissipationModel = new DissipationModelPow(dissipationModelParameters[0]);
                } else {
                    std::cerr << "[ERROR] dissipation model parameters for power dissipation must be one: aborting"<<std::endl;
                    return 1;
                }
            } else {
                std::cerr << "[ERROR] dissipation model parameters for power dissipation was not set: setting to default (2)"<<std::endl;
                dissipationModel = new DissipationModelPow(2);
            }
        } else if(dissipationModelName == "random"){
            if (vm.count("dissipationModelParameters")) {
                std::cout << "[LOG] dissipation model parameters were declared to be "
            << vm["dissipationModelParameters"].as<std::vector<double>>()[0] << " & " << vm["dissipationModelParameters"].as<std::vector<double>>()[1] << ".\n";
                std::vector<double> dissipationModelParameters = vm["dissipationModelParameters"].as<std::vector<double>>();
                if(dissipationModelParameters.size() == 2){
                    dissipationModel = new DissipationModelRandom(dissipationModelParameters[0],dissipationModelParameters[1]);
                } else {
                    std::cerr << "[ERROR] dissipation model parameters for random dissipation must be two: aborting"<<std::endl;
                    return 1;
                }
            } else {
                std::cerr << "[ERROR] dissipation model parameters for random dissipation was not set: aborting"<<std::endl;
                return 1;
            }
        } else if(dissipationModelName == "scaled"){
            if (vm.count("dissipationModelParameters")) {
                std::cout << "[LOG] dissipation model parameters were declared to be "
            << vm["dissipationModelParameters"].as<std::vector<double>>()[0] << ".\n";
                std::vector<double> dissipationModelParameters = vm["dissipationModelParameters"].as<std::vector<double>>();
                if(dissipationModelParameters.size() == 1){
                    dissipationModel = new DissipationModelScaled([dissipationModelParameters](double time)->double{return dissipationModelParameters[0];});
                } else {
                    std::cerr << "[ERROR] dissipation model parameters for scaled dissipation must be one: aborting"<<std::endl;
                    return 1;
                }
            } else {
                std::cerr << "[ERROR] dissipation model parameters for scaled dissipation was not set: setting to default 0.5 costant"<<std::endl;
                dissipationModel = new DissipationModelScaled();
            }
        } else if(dissipationModelName == "periodic"){
            // std::cout << "Not yet implemented\n";
            // return 1;
            // if (vm.count("dissipationModelParameters")) {
            //     std::vector<double> dissipationModelParameters = vm["dissipationModelParameters"].as<std::vector<double>>();
            //     if (dissipationModelParameters.size() == 3) {
            //         std::cout << "[LOG] dissipation model parameters were set to "
            //         << dissipationModelParameters[0] << " & " << dissipationModelParameters[1] << " & " << dissipationModelParameters[2] << ".\n";
            //         dissipationModel = new DissipationModelPeriodic(dissipationModelParameters[0],dissipationModelParameters[1],dissipationModelParameters[2]);
            //     } else {
            //         std::cerr << "[ERROR] dissipation model parameters for periodic dissipation must be three: aborting"<<std::endl;
            //         return 1;
            //     }
                
                
            // }
            if (vm.count("dissipationModelParameters")) {
                std::vector<double> dissipationModelParameters = vm["dissipationModelParameters"].as<std::vector<double>>();
                if (dissipationModelParameters.size() == 3) {
                    std::cout << "[LOG] dissipation model parameters were set to Amplitude:"
                    << dissipationModelParameters[0] << " & period:" << dissipationModelParameters[1] << " & phase: " << dissipationModelParameters[2] << ".\n";
                    dissipationModel = new DissipationModelScaled([dissipationModelParameters](double time)->double{return dissipationModelParameters[0]*sin(2*arma::datum::pi/dissipationModelParameters[1]*time + dissipationModelParameters[2]);});
                } else {
                    std::cerr << "[ERROR] dissipation model parameters for periodic dissipation must be three for amplitude, period and phase: aborting"<<std::endl;
                    return 1;
                }
                
                
            } else {
                std::cerr << "[ERROR] dissipation model parameters for periodic dissipation was not set: aborting"<<std::endl;
                return 1;
            }
        } 
    } else { //dissipation model set to default (none)
        std::cout << "[LOG] dissipation model was not set. set to default (none)\n";
        dissipationModel = new DissipationModel();
    }
    //end program options section

    std::map<std::string,std::string> ensembletoEntrez = getEnsembletoEntrezidMap();
    if(ensembleGeneNames){
        std::cout <<"[LOG] mapping ensemble gene names to entrez ids"<<std::endl;
    }
    auto namesAndEdges = edgesFileToEdgesListAndNodesByName(filename);
    std::vector<std::string> metapathwayNodes = namesAndEdges.first;
    WeightedEdgeGraph *metapathway = new WeightedEdgeGraph(metapathwayNodes);
    for(auto edge = namesAndEdges.second.cbegin() ; edge != namesAndEdges.second.cend(); edge++ ){
        metapathway->addEdge(std::get<0> (*edge), std::get<1> (*edge) ,std::get<2>(*edge) );
    }


    auto logFolds = logFoldChangeMatrixToCellVectors(cellLogFoldMatrixFilename,metapathwayNodes,ensembleGeneNames);
    std::vector<std::string> geneslogfoldNames = std::get<0>(logFolds);
    std::vector<std::string> cellTypes = std::get<1>(logFolds);
    Computation** cellComputations = new Computation*[cellTypes.size()];
    for(uint i = 0; i < cellTypes.size();i++){
        std::vector<double> inputCelllogfold = std::get<2>(logFolds)[i];
        Computation* tmpCompPointer = new Computation(cellTypes[i],inputCelllogfold,metapathway,metapathwayNodes);  //TODO order the genes directly or use the names and set them one by one 
        tmpCompPointer->setDissipationModel(dissipationModel);
        cellComputations[i] = tmpCompPointer;
        //No inverse computation with the augmented pathway since virtual nodes edges are not yet inserted
        cellComputations[i]->augmentMetapathwayNoComputeInverse(cellTypes);
    }
    auto allFilesInteraction = get_all(celltypesInteractionFoldername,".tsv");
    for(auto cellInteractionFilename = allFilesInteraction.cbegin() ; cellInteractionFilename != allFilesInteraction.cend() ; cellInteractionFilename++){
        auto cellInteractionsEdges = cellInteractionFileToEdgesListAndNodesByName(*cellInteractionFilename,ensembleGeneNames);
        //TODO insert edges to the correspondent cell metapathway
        #pragma omp parallel for
        for (uint i = 0; i < cellTypes.size();i++) {
            if(cellInteractionsEdges.contains(cellTypes[i])){
                cellComputations[i]->addEdges(cellInteractionsEdges[cellTypes[i]]);
                //cellComputations[i]->freeAugmentedGraphs();  //REMOVE
            }
        }
    }


    //freeing some data structures inside computation to consume less RAM
    std::vector<std::vector<std::string>> cellToNodeNames = std::vector<std::vector<std::string>>(cellTypes.size(),std::vector<std::string>());
    for(uint i = 0; i < cellTypes.size();i++ ){
        cellToNodeNames[i] = cellComputations[i]->getAugmentedMetapathway()->getNodeNames();  
        cellComputations[i]->freeAugmentedGraphs();  
    }

    // EndCelltype -> (sourceCellType -> value)

    uint iterationIntercell = 0;
    while(iterationIntercell++ < intercellIterations){
        //computation of perturbation
        uint iterationIntracell = 0;
        //intracell iteration with no passing of values to the virtual nodes
        while (iterationIntracell++ < intracellIterations) {
            #pragma omp parallel for
            for(uint i = 0; i < cellTypes.size(); i++){
                std::vector<std::string> nodeNames = cellToNodeNames[i];
                std::cout << "[LOG] computation of perturbation for iteration ("+ std::to_string(iterationIntercell) + ") for cell (" + cellTypes[i]<<std::endl; 
                //std::vector<double> outputValues = cellComputations[i]->computeAugmentedPerturbation();
                //std::vector<double> outputValues = cellComputations[i]->computeAugmentedPerturbationSaturated();
                std::vector<double> outputValues = cellComputations[i]->computeAugmentedPerturbationDissipatedBeforeCompute(iterationIntracell); // TODO check if iteration intracell should be multiplied by iteration intercell
            }
            //save output values
            for(uint i = 0; i < cellTypes.size(); i++){
                std::vector<std::string> nodeNames = cellToNodeNames[i];
                //TODO change how to save files to get more information about intracell and intercell iterations
                saveNodeValues(outputFoldername, iterationIntercell*intracellIterations + iterationIntracell, cellTypes[i], cellComputations[i]->getOutputAugmented(), nodeNames,ensembleGeneNames);
            }
        }
        //update input
        for(uint i = 0; i < cellTypes.size(); i++){
            std::cout << "[LOG] update input for iteration ("+ std::to_string(iterationIntercell) + ") for cell (" + cellTypes[i]<<std::endl;
            cellComputations[i]->updateInput(std::vector<double>(),true);
        }
        //update input with virtual node values update

        for (uint i = 0; i < cellTypes.size(); i++) {
            //queuesCellTypes[i] = cellComputations[i]->computeAugmentedPerturbation();
            //TODO when computation will be done in parallel, this step should wait for all the computations of the other adjacent cells to finish
            for(uint j = 0; j < cellTypes.size(); j++){
                if(i==j){
                    if(sameCellCommunication) cellComputations[i]->setInputVinForCell(cellTypes[j], cellComputations[j]->getVirtualOutputForCell(cellTypes[i]));
                } else {
                    cellComputations[i]->setInputVinForCell(cellTypes[j], cellComputations[j]->getVirtualOutputForCell(cellTypes[i]));
                }
            }
        }
        
    }
    


    //cleaning memory
    //delete metapathway
    

    
    return 0;
}