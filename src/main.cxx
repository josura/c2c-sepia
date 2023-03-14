#include <boost/program_options/value_semantic.hpp>
#include <iostream>
#include <boost/program_options.hpp>
#include <tuple>
#include <vector>
#include "Computation.h"
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
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    std::string filename,celltypesFilename,celltypesInteractionFoldername,cellLogFoldMatrixFilename;

    if (vm.count("help")) {
        //printHelp();
        std::cout << desc << std::endl;
        return 1;
    }

    if (vm.count("fMETAPATHWAY")) {
        std::cout << "[LOG] file for the metapathway was set to " 
    << vm["fMETAPATHWAY"].as<std::string>() << ".\n";
        filename = vm["fMETAPATHWAY"].as<std::string>();
    } else {
        std::cout << "[ERROR] fMETAPATHWAY file was not set. Aborting\n";
        return 1;
    }
    if (vm.count("fLogfoldPerCell")) {
        std::cout << "[LOG] file for the logfoldPerCell matrix was set to " 
    << vm["fLogfoldPerCell"].as<std::string>() << ".\n";
        cellLogFoldMatrixFilename = vm["fLogfoldPerCell"].as<std::string>();
    } else {
        std::cout << "[ERROR] fLogfoldPerCell file was not set. Aborting\n";
        return 1;
    }
    if (vm.count("dirCellInteraction")) {
        std::cout << "[LOG] folder for the cell interactions was set to " 
    << vm["dirCellInteraction"].as<std::string>() << ".\n";
        celltypesInteractionFoldername = vm["dirCellInteraction"].as<std::string>();
    } else {
        std::cout << "[LOG] dirCellInteraction folder was not set. computing without taking into account cell interactions\n";
        //TODO
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
    std::vector<Computation> cellComputations;
    for(uint i = 0; i < cellTypes.size();i++){
        std::vector<double> inputCelllogfold = std::get<2>(logFolds)[i];
        cellComputations.push_back(Computation(cellTypes[i],inputCelllogfold,metapathway,metapathwayNodes));  //TODO order the genes directly or use the names and set them one by one
        //TODO augment the metapathway, I am scared since I do not have a lot of memory
        cellComputations[i].augmentMetapathway(cellTypes);
    }
    auto allFilesInteraction = get_all(celltypesInteractionFoldername,".tsv");
    for(auto cellInteractionFilename = allFilesInteraction.cbegin() ; cellInteractionFilename != allFilesInteraction.cend() ; cellInteractionFilename++){
        auto cellInteractionsEdges = cellInteractionFileToEdgesListAndNodesByName(*cellInteractionFilename,ensembleGeneNames);
        //TODO insert edges to the correspondent cell metapathway
        for (uint i = 0; i < cellTypes.size();i++) {
            if(cellInteractionsEdges.contains(cellTypes[i])){
                cellComputations[i].addEdges(cellInteractionsEdges[cellTypes[i]]);
                cellComputations[i].freeAugmentedGraphs();  //REMOVE
            }
        }
    }

    //freeing some data structures inside computation to consume less RAM
    for(uint i = 0; i < cellComputations.size();i++ ){
        cellComputations[i].freeAugmentedGraphs();
    }


    //cleaning memory
    //delete metapathway
    

    
    return 0;
}