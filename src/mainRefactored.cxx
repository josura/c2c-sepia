#include <boost/program_options/value_semantic.hpp>
#include <iostream>
#include <boost/program_options.hpp>
#include <map>
#include <ostream>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <vector>
#include "Computation.h"
#include "ConservationModel.h"
#include "DissipationModel.h"
#include "DissipationModelPow.h"
#include "DissipationModelRandom.h"
#include "DissipationModelScaled.h"
#include "WeightedEdgeGraph.h"
#include "utilities.h"
#include "CustomScalingFunctions.h"


void printHelp(){
    //TODO fix this help
    std::cout << "usage: ./c2c-sepia --fUniqueGraph <graph>.tsv --fInitialPerturbationPerType <initialPerturbationPerType>.tsv [<subtypes>.txt] --typeInteractionFolder <TypeInteractionFolder>(containing .tsv files)]\nFILE STRUCTURE SCHEMA:\ngraph.tsv\nstart\tend\tweight\n<gene1>\t<gene2>\t <0.something>\n...\n\n\ninitialPerturbationPerType.tsv\n\ttype1\ttype2\t...\ttypeN\ngene1\t<lfc_type1:gene1>\t<lfc_type2:gene1>\t...\t<lfc_typeN:gene1>\ngene1\t<lfc_type1:gene2>\t<lfc_type2:gene2>\t...\t<lfc_typeN:gene2>\n...\n\n\ntypesInteraction.tsv\nstartType:geneLigand\tendType:geneReceptor\tweight\n<type1:geneLigand>\t<type2:genereceptor>\t <0.something>\n...\n\n\nsubtypes.txt\ntype1\ntype3\n..."<<std::endl;
    std::cout << "LEGEND:\n <> := placeholder for the name of the file\n[] := optional\n{} := at least one"<<std::endl;
}

int main(int argc, char** argv ) {
    //program options
    bool ensembleGeneNames=false;
    bool sameTypeCommunication=false;
    bool saturation=false;
    bool conservateInitialNorm=false;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    //TODO implement subtypes
    desc.add_options()
        ("help", "() print help section")//<initialPerturbationPerType>.tsv [<subtypes>.txt] [<typesInteraction>.tsv]\nFILE STRUCTURE SCHEMA:\ngraph.tsv\nstart end weight\n<gene1> <gene2>  <0.something>\n...\n\n\ninitialPerturbationPerType.tsv\n type1 type2 ... typeN\ngene1 <lfc_type1:gene1> <lfc_type2:gene1> ... <lfc_typeN:gene1>\ngene1 <lfc_type1:gene2> <lfc_type2:gene2> ... <lfc_typeN:gene2>\n...\n\n\ntypesInteraction.tsv\nstartType:geneLigand endType:geneReceptor weight\n<type1:geneLigand> <type2:genereceptor>  <0.something>\n...\n\n\nsubtypes.txt\ntype1\ntype3\n...")
        ("fUniqueGraph", po::value<std::string>(), "(string) graph filename, for an example graph see in resources. NOTE: if this option is chosen, graphsFilesFolder cannot be used")
        ("fInitialPerturbationPerType", po::value<std::string>(), "(string) initialPerturbationPerType matrix filename, for an example see in data")
        ("subtypes", po::value<std::string>(), "subtypes filename, for an example see in data TODO")
        ("initialPerturbationPerTypeFolder", po::value<std::string>(), "(string) initialPerturbationPerType folder, for an example see in data TODO")
        ("typeInteractionFolder", po::value<std::string>(), "(string) directory for the type interactions, for an example see in data")
        ("ensembleGeneNames",po::bool_switch(&ensembleGeneNames),"() use ensemble gene names, since the graph used in resources have entrez_ids, a map will be done from ensemble to entrez, the map is available in resources")
        ("sameTypeCommunication",po::bool_switch(&sameTypeCommunication),"() use same type communication, since it is not permitted as the standard definition of the model, this adds a virtual node for the same type type")
        ("outputFolder",po::value<std::string>()->required(),"(string) output folder for output of the algorithm at each iteration")
        ("intertypeIterations",po::value<uint>(),"(positive integer) number of iterations for intertype communication")
        ("intratypeIterations",po::value<uint>(),"(positive integer) number of iterations for intratype communication")
        ("timestep",po::value<double>(),"timestep to use for the iteration, the final time is iterationIntracell*iterationIntercell*timestep")
        ("dissipationModel",po::value<std::string>(),"(string) the dissipation model for the computation, available models are: 'none (default)','power','random','periodic','scaled' and 'custom'")
        ("dissipationModelParameters",po::value<std::vector<double>>()->multitoken(),"(string) the parameters for the dissipation model, for the power dissipation indicate the base, for the random dissipation indicate the min and max value, for the periodic dissipation indicate the period")
        ("graphsFilesFolder",po::value<std::string>(),"(string) graphs (pathways or other types of graphs) file folder TODO implement different graphs loading")
        ("conservationModel",po::value<std::string>(),"(string) the conservation model used for the computation, available models are: 'none (default)','scaled','random' and 'custom' ")
        ("conservationModelParameters", po::value<std::vector<double>>()->multitoken(),"(vector<double>) the parameters for the dissipation model, for the scaled parameter the constant used to scale the conservation final results, in the case of random the upper and lower limit (between 0 and 1)")
        ("saturation",po::bool_switch(&saturation),"use saturation of values, default to 1, if another value is needed, use the saturationTerm")
        ("saturationTerm",po::value<double>(),"defines the limits of the saturation [-saturationTerm,saturationTerm]")
        ("conservateInitialNorm",po::bool_switch(&conservateInitialNorm), "conservate the initial euclidean norm of the perturbation values, that is ||Pn|| <= ||Initial||, default to false")
    ;
    //TODO add additional boolean parameter to control if the graph names are not genes and the algorithm should use the graph names directly, no conversion or mapping

    

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    std::string filename,subtypesFilename,typesFilename,typesInteractionFoldername,typesInitialPerturbationMatrixFilename,graphsFilesFolder, typeInitialPerturbationFolderFilename,outputFoldername;
    uint intertypeIterations,intratypeIterations;
    DissipationModel* dissipationModel = nullptr;
    ConservationModel* conservationModel = nullptr;
    double timestep = 1;

    if (vm.count("help")) {
        //printHelp();
        std::cout << desc << std::endl;
        return 1;
    }


    //controls over impossible configurations
    if(vm.count("fUniqueGraph") == 0 && vm.count("graphsFilesFolder") == 0){
        //no unique graph of folder of the graphs was set
        std::cout << "[ERROR] no unique graph filename or folder was set to get the graphs, set one "<<std::endl;
        return 1;
    }

    if(vm.count("fInitialPerturbationPerType") == 0 && vm.count("initialPerturbationPerTypeFolder") == 0){
        //no way of getting the initial perturbation values
        std::cout << "[ERROR] no matrix for the initial values was passed as filename or single vector in files contained in the folder specified was set, set one "<<std::endl;
        return 1;
    }

    if(vm.count("fInitialPerturbationPerType") && vm.count("graphsFilesFolder")){
        //unstable configuration of different graphs and single matrix with the same nodes
        std::cout << "[WARNING] unstable configuration of different graphs and a single matrix with the initial perturbations"<<std::endl;
    }

    if(saturation && conservateInitialNorm){
        std::cerr << "[ERROR] saturation and conservateInitialNorm cannot be both true, aborting"<<std::endl;
        return 1;
    }


    if(vm.count("graphsFilesFolder") && vm.count("fUniqueGraph")){
        std::cout << "[ERROR] fUniqueGraph and graphsFilesFolder were both set. Aborting\n";
        return 1;
    }
    if(vm.count("initialPerturbationPerTypeFolder") && vm.count("fInitialPerturbationPerType")){
        std::cout << "[ERROR] fInitialPerturbationPerType and initialPerturbationPerTypeFolder were both set. Aborting\n";
        return 1;
    }

    // reading the parameters

    if (vm.count("intertypeIterations")) {
        std::cout << "[LOG] iterations intertype set to " 
    << vm["intertypeIterations"].as<std::string>() << ".\n";
        intertypeIterations = vm["intertypeIterations"].as<uint>();
    } else {
        std::cout << "[LOG] iterations intertype not set, set to default: 10 iterations \n";
        intertypeIterations = 10;
    }

    if (vm.count("intratypeIterations")) {
        std::cout << "[LOG] iterations intratype set to " 
    << vm["intratypeIterations"].as<std::string>() << ".\n";
        intratypeIterations = vm["intratypeIterations"].as<uint>();
    } else {
        std::cout << "[LOG] iterations intratype not set, set to default: 5 iterations \n";
        intratypeIterations = 5;
    }


    if(vm.count("timestep")){
        std::cout << "[LOG] timestep set to " 
    << vm["timestep"].as<std::string>() << ".\n";
        timestep = vm["timestep"].as<double>();
    } else {
        std::cout << "[LOG] timestep not set, set to default (1)"<<std::endl;
    }


    if (vm.count("fUniqueGraph")) {
        std::cout << "[LOG] file for the graph was set to " 
    << vm["fUniqueGraph"].as<std::string>() << ".\n";
        filename = vm["fUniqueGraph"].as<std::string>();
        if(!fileExistsPath(filename)){
            std::cerr << "[ERROR] file for the graph do not exist: aborting"<<std::endl;
            return 1;
        }
    } else if(vm.count("graphsFilesFolder")){
        std::cout << "[LOG] folder for the graphs was set to " 
    << vm["graphsFilesFolder"].as<std::string>() << ".\n";
        graphsFilesFolder = vm["graphsFilesFolder"].as<std::string>();
        if(!folderExists(graphsFilesFolder)){
            std::cerr << "[ERROR] folder for the graphs do not exist: aborting"<<std::endl;
            return 1;
        }
    }
    if (vm.count("fInitialPerturbationPerType")) {
        std::cout << "[LOG] file for the initialPerturbationPerType matrix was set to " 
    << vm["fInitialPerturbationPerType"].as<std::string>() << ".\n";
        typesInitialPerturbationMatrixFilename = vm["fInitialPerturbationPerType"].as<std::string>();
        if(!fileExistsPath(typesInitialPerturbationMatrixFilename)){
            std::cerr << "[ERROR] file for the initialPerturbationPerType does not exist: aborting"<<std::endl;
            return 1;
        }
    } else if (vm.count("initialPerturbationPerTypeFolder")) {
        std::cout << "[LOG] folder for the initialPerturbationPerType was set to "
    << vm["initialPerturbationPerTypeFolder"].as<std::string>() << ".\n";
        typeInitialPerturbationFolderFilename = vm["initialPerturbationPerTypeFolder"].as<std::string>();
        if(!folderExists(typeInitialPerturbationFolderFilename)){
            std::cerr << "[ERROR] folder for the initialPerturbationPerType do not exist: aborting"<<std::endl;
            return 1;
        }
    }

    if (vm.count("typeInteractionFolder")) {
        std::cout << "[LOG] folder for the type interactions was set to " 
    << vm["typeInteractionFolder"].as<std::string>() << ".\n";
        typesInteractionFoldername = vm["typeInteractionFolder"].as<std::string>();
        if(!folderExists(typesInteractionFoldername)){
            std::cerr << "[ERROR] folder for the type interactions do not exist: aborting"<<std::endl;
            return 1;
        }
    } else {
        std::cout << "[LOG] typeInteractionFolder folder was not set. computing without taking into account type interactions\n";
        //TODO
    }
    if (vm.count("outputFolder")) {
        std::cout << "[LOG] output folder  was set to " 
    << vm["outputFolder"].as<std::string>() << ".\n";
        outputFoldername = vm["outputFolder"].as<std::string>();
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
            dissipationModel = new DissipationModelScaled([](double time)->double{return 0;});
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
        } else if(dissipationModelName == "custom"){
            //control if custom function for dissipation returns double and takes a single parameter as double
            std::cout << "[LOG] dissipation model was set to custom, HIGH RISK OF FAILURE\n " << std::endl;
            dissipationModel = new DissipationModelScaled(getDissipationScalingFunction());
        } else {
            std::cerr << "[ERROR] dissipation model scale function is not any of the types. Conservation model scale functions available are none(default), scaled, random and custom \n";
            return 1;
        }
    } else { //dissipation model set to default (none)
        std::cout << "[LOG] dissipation model was not set. set to default (none)\n";
        dissipationModel = new DissipationModelScaled([](double time)->double{return 0;});
    }


    if (vm.count("conservationModel")) {
        std::cout << "[LOG] conservation model was set to "
    << vm["conservationModel"].as<std::string>() << ".\n";
        std::string conservationModelName = vm["conservationModel"].as<std::string>();
        if(conservationModelName == "none"){
            std::cout << "[LOG] conservation model set to default (none)\n";
            conservationModel = new ConservationModel([](double time)->double{return 0;});
        } else if (conservationModelName == "scaled"){
            if (vm.count("conservationModelParameters")) {
                std::cout << "[LOG] conservation model parameters were declared to be "
            << vm["conservationModelParameters"].as<std::vector<double>>()[0] << ".\n";
                std::vector<double> conservationModelParameters = vm["conservationModelParameters"].as<std::vector<double>>();
                if(conservationModelParameters.size() == 1){
                    conservationModel = new ConservationModel([conservationModelParameters](double time)->double{return conservationModelParameters[0];});
                } else {
                    std::cerr << "[ERROR] conservation model parameters for scaled conservation must be one parameter: aborting"<<std::endl;
                    return 1;
                }
            } else {
                std::cerr << "[ERROR] conservation model parameters for scaled conservation was not set: setting to default 0.5 costant"<<std::endl;
                conservationModel = new ConservationModel();
            }
        } else if (conservationModelName == "random"){
            if (vm.count("conservationModelParameters")) {
                std::cout << "[LOG] conservation model parameters were declared to be "
            << vm["conservationModelParameters"].as<std::vector<double>>()[0] << " & " << vm["conservationModelParameters"].as<std::vector<double>>()[1] << ".\n";
                std::vector<double> conservationModelParameters = vm["conservationModelParameters"].as<std::vector<double>>();
                if(conservationModelParameters.size() == 2){
                    //control if lower and upper limits of the random values are within 0 and 1
                    if( (conservationModelParameters[0] < 0) || (conservationModelParameters[0] > 1) || (conservationModelParameters[1] < 0) || (conservationModelParameters[1] > 1) || (conservationModelParameters[0] > conservationModelParameters[1]) ){
                        std::cerr << "[ERROR] conservation model parameters for random conservation must be between 0 and 1 and must be a < b: aborting"<<std::endl;
                        return 1;
                    }
                    conservationModel = new ConservationModel([conservationModelParameters](double time)->double{return randomRealNumber(conservationModelParameters[0],conservationModelParameters[1]);});
                } else {
                    std::cerr << "[ERROR] conservation model parameters for random conservation must be two: aborting"<<std::endl;
                    return 1;
                }
            } else {
                std::cerr << "[ERROR] conservation model parameters for random conservation was not set: aborting"<<std::endl;
                return 1;
            }
        } else if(conservationModelName == "custom"){
            //control if custom function for conservation returns double and takes a single parameter as double
            std::cout << "[LOG] conservation model was set to custom, HIGH RISK OF FAILURE\n " << std::endl;
            conservationModel = new ConservationModel(getConservationScalingFunction());
        } else {
            std::cerr << "[ERROR] conservation model scale function is not any of the types. Conservation model scale functions available are none(default), scaled, random and custom \n";
            return 1;
        }
    } else {
        std::cout << "[LOG] conservation model was not set. set to default (none)\n";
        conservationModel = new ConservationModel([](double time)->double{return 0;});
    }

    //logging if saturation is set and saturation parameters are set
    if (saturation) {
        if(vm.count("saturationTerm") == 0){
            std::cout << "[LOG] saturation term not specified, using the interval [-1,1]"<<std::endl;
        } else if(vm.count("saturationTerm") == 1){
            double saturationTerm = vm["saturationTerm"].as<double>();
            std::cout << "[LOG] saturation term specified, using the interval [-" << saturationTerm << "," << saturationTerm << "]"<<std::endl;
        } else {
            std::cerr << "[ERROR] saturation term specified more than once, possibility of using more values not yet implemented: aborting"<<std::endl;
            return 1;
        }
    }
    //end program options section

    std::map<std::string,std::string> ensembletoEntrez = getEnsembletoEntrezidMap();
    if(ensembleGeneNames){
        std::cout <<"[LOG] mapping ensemble gene names to entrez ids"<<std::endl;
    }
    //take the types before with another function TODO define function
    std::vector<std::string> types;
    if(vm.count("fUniqueGraph")){
        if(vm.count("fInitialPerturbationPerType")){
            types = getTypesFromMatrixFile(typesInitialPerturbationMatrixFilename);

        } else if (vm.count("initialPerturbationPerTypeFolder")){
            types = getTypesFromFolderFileNames(typeInitialPerturbationFolderFilename);
        } else {
            std::cerr << "[ERROR] no initial perturbation file or folder specified: aborting"<<std::endl;
            return 1;
        }
    } else if (vm.count("graphsFilesFolder")) {
        types = getTypesFromFolderFileNames(graphsFilesFolder);
    } else {
        std::cerr << "[ERROR] no graph file or folder specified: aborting"<<std::endl;
        return 1;
    }

    std::vector<std::string> subtypes;
    if(vm.count("subtypes")){
        std::cout << "[LOG] subtypes filename set to "
    << vm["subtypes"].as<std::string>() << ".\n";
        subtypesFilename = vm["subtypes"].as<std::string>();
        subtypes = getVectorFromFile<std::string>(subtypesFilename);
    }else{
        std::cout << "[LOG] subtypes filename not set, set to default: all types \n";
        subtypes = types;
    }
    
    //filter types with the subtypes
    std::vector<std::string> typesFiltered = vectorsIntersection(types, subtypes);
    if (typesFiltered.size() == 0) {
        std::cerr << "[ERROR] no types in common between the types and subtypes: aborting"<<std::endl;
        return 1;
    }
    
    
    //use the number of types to allocate an array of pointers to contain the graph for every type
    WeightedEdgeGraph **graphs = new WeightedEdgeGraph*[types.size()];
    std::vector<std::vector<std::string>> graphsNodes;
    std::vector<std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>>> namesAndEdges;
    if(vm.count("fUniqueGraph")){
        namesAndEdges.push_back(edgesFileToEdgesListAndNodesByName(filename));
        graphsNodes.push_back(namesAndEdges[0].first);
        graphs[0] = new WeightedEdgeGraph(graphsNodes[0]);
        for(uint i = 1; i < types.size(); i++){
            namesAndEdges.push_back(namesAndEdges[0]);
            graphsNodes.push_back(namesAndEdges[0].first);
            graphs[i] = graphs[0];
        }
    } else if (vm.count("graphsFilesFolder")) {
        // TODO get the nodes from the single files
        auto allGraphs = edgesFileToEdgesListAndNodesByNameFromFolder(graphsFilesFolder);
        auto typesFromFolder = allGraphs.first;
        if(typesFromFolder.size() != types.size()){
            std::cerr << "[ERROR] types from folder and types from file do not match: aborting"<<std::endl;
            return 1;
        }
        for (uint i = 0; i<typesFromFolder.size(); i++){
            if(typesFromFolder[i] != types[i]){
                std::cerr << "[ERROR] types from folder and types from file do not match: aborting"<<std::endl;
                return 1;
            }
        }
        namesAndEdges = allGraphs.second;
        for(uint i = 0; i < types.size(); i++){
            graphsNodes.push_back(namesAndEdges[i].first);
            graphs[i] = new WeightedEdgeGraph(graphsNodes[i]);
        }
    } 

    //add the edges to the graphs
    if(vm.count("fUniqueGraph")){
        for(auto edge = namesAndEdges[0].second.cbegin() ; edge != namesAndEdges[0].second.cend(); edge++ ){
            graphs[0]->addEdge(std::get<0> (*edge), std::get<1> (*edge) ,std::get<2>(*edge) );
        }
    } else if (vm.count("graphsFilesFolder")) {
        for(uint i = 0; i < types.size(); i++){
            for(auto edge = namesAndEdges[i].second.cbegin() ; edge != namesAndEdges[i].second.cend(); edge++ ){
                graphs[i]->addEdge(std::get<0> (*edge), std::get<1> (*edge) ,std::get<2>(*edge) );
            }
        }
    }


    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::vector<double>>> initialValues;
    std::vector<std::vector<double>> inputInitials;
    if(vm.count("fInitialPerturbationPerType")){
        std::cout << "[LOG] initial perturbation per type specified, using the file "<<typesInitialPerturbationMatrixFilename<<std::endl;
        if(subtypes.size()==0){
            std::cout << "[LOG] no subcelltypes specified, using all the celltypes in the log fold matrix"<<std::endl;
            initialValues = logFoldChangeMatrixToCellVectors(typesInitialPerturbationMatrixFilename,graphsNodes[0],ensembleGeneNames);
        } else {
            std::cout << "[LOG] subcelltypes specified, using only the celltypes in the log fold matrix that are in the list"<<std::endl;
            initialValues = logFoldChangeMatrixToCellVectors(typesInitialPerturbationMatrixFilename,graphsNodes[0],subtypes,ensembleGeneNames);
        }
    } else if (vm.count("initialPerturbationPerTypeFolder")){
        std::cout << "[LOG] initial perturbation per type specified, using the folder "<<typeInitialPerturbationFolderFilename<<std::endl;
        initialValues = logFoldChangeCellVectorsFromFolder(typeInitialPerturbationFolderFilename,graphsNodes,subtypes,ensembleGeneNames);
    } else {
        std::cerr << "[ERROR] no initial perturbation file or folder specified: aborting"<<std::endl;
        return 1;
    }
    std::vector<std::string> initialNames = std::get<0>(initialValues);
    inputInitials = std::get<2>(initialValues);
    std::vector<std::string> typesFromValues = std::get<1>(initialValues);
    //TODO understand if types from values should be the same as the types from the graphs since values could be specified for a subset of the types
    //this condition should take into account the intersection of the types and the subtypes
    if(typesFromValues.size() == 0){
        std::cerr << "[ERROR] types from the initial values folder are 0, control if the types are the same to the one specified in the matrix, in the graphs folder and in the subtypes: aborting"<<std::endl;
        std::cerr << "[ERROR] types specified(subtypes): ";
        for(auto type: subtypes)
            std::cerr << type << " ";
        std::cerr << std::endl;
        std::cerr << "[ERROR] types from file(from graphs folder or from matrix): ";
        for(auto type: types)
            std::cerr << type << " ";
        std::cerr << std::endl;
        std::cerr << "[ERROR] types from values(from initial values folder or from values matrix) intersected with subtypes: ";
        for(auto type: typesFromValues)
            std::cerr << type << " ";
        std::cerr << std::endl;
        return 1;
    }
    auto indexMapGraphTypesToValuesTypes = get_indexmap_vector_values_full(types, typesFromValues);
    if(indexMapGraphTypesToValuesTypes.size() == 0){
        std::cerr << "[ERROR] types from folder and types from file do not match even on one instance: aborting"<<std::endl;
        return 1;
    }

    Computation** typeComputations = new Computation*[typesFiltered.size()];
    int indexComputation = 0;
    std::vector<int> typesIndexes = std::vector<int>(types.size(),-1); 
    std::vector<int> invertedTypesIndexes = std::vector<int>(typesFiltered.size(),-1); 
    for(uint i = 0; i < types.size();i++){
        if(vectorContains(typesFiltered, types[i])){
            if(indexMapGraphTypesToValuesTypes[i] == -1){
                std::cout << "[LOG] type "<<types[i]<<" not found in the initial perturbation files, using zero vector as input"<<std::endl;
                std::vector<double> input = std::vector<double>(graphsNodes[i].size(),0);
                Computation* tmpCompPointer = new Computation(types[i],input,graphs[i],graphsNodes[i]);  //TODO order the genes directly or use the names and set them one by one 
                tmpCompPointer->setDissipationModel(dissipationModel);
                tmpCompPointer->setConservationModel(conservationModel);
                typeComputations[indexComputation] = tmpCompPointer;
                //No inverse computation with the augmented pathway since virtual nodes edges are not yet inserted
                typeComputations[indexComputation]->augmentMetapathwayNoComputeInverse(typesFiltered);
            } else {
                int index = indexMapGraphTypesToValuesTypes[i];
                std::vector<double> input = inputInitials[index];
                Computation* tmpCompPointer = new Computation(types[i],input,graphs[i],graphsNodes[i]);  //TODO order the genes directly or use the names and set them one by one 
                tmpCompPointer->setDissipationModel(dissipationModel);
                tmpCompPointer->setConservationModel(conservationModel);
                typeComputations[indexComputation] = tmpCompPointer;
                //No inverse computation with the augmented pathway since virtual nodes edges are not yet inserted
                typeComputations[indexComputation]->augmentMetapathwayNoComputeInverse(typesFiltered);
            }
            typesIndexes[i] = indexComputation;
            invertedTypesIndexes[indexComputation] = i;
            indexComputation++;
        }
    }
    std::vector<std::vector<std::string>> typeToNodeNames = std::vector<std::vector<std::string>>(typesFiltered.size(),std::vector<std::string>());
    for(uint i = 0; i < typesFiltered.size();i++ ){
        typeToNodeNames[i] = typeComputations[i]->getAugmentedMetapathway()->getNodeNames();    
    }
    auto allFilesInteraction = get_all(typesInteractionFoldername,".tsv");
    for(auto typeInteractionFilename = allFilesInteraction.cbegin() ; typeInteractionFilename != allFilesInteraction.cend() ; typeInteractionFilename++){
        std::map<std::string, std::vector<std::tuple<std::string, std::string, double>>> typeInteractionsEdges;
        if (subtypes.size() == 0) {
            typeInteractionsEdges  = cellInteractionFileToEdgesListAndNodesByName(*typeInteractionFilename,ensembleGeneNames);
        } else {
            typeInteractionsEdges = cellInteractionFileToEdgesListAndNodesByName(*typeInteractionFilename, subtypes, ensembleGeneNames);
        }
        //TODO insert edges to the correspondent type graph
        #pragma omp parallel for
        for (uint i = 0; i < types.size();i++) {
            if(typeInteractionsEdges.contains(types[i]) && typesIndexes[i] != -1){
                typeComputations[typesIndexes[i]]->addEdges(typeInteractionsEdges[types[i]]);
                //typeComputations[i]->freeAugmentedGraphs();
            }
        }
    }


    //freeing some data structures inside computation to consume less RAM
    // std::vector<std::vector<std::string>> typeToNodeNames = std::vector<std::vector<std::string>>(types.size(),std::vector<std::string>());
    // for(uint i = 0; i < types.size();i++ ){
    //     typeToNodeNames[i] = typeComputations[i]->getAugmentedMetapathway()->getNodeNames();  
    //     typeComputations[i]->freeAugmentedGraphs();  
    // }

    // EndTypetype -> (sourceTypeType -> value)

    uint iterationIntertype = 0;
    while(iterationIntertype < intertypeIterations){
        //computation of perturbation
        uint iterationIntratype = 0;
        //intratype iteration with no passing of values to the virtual nodes
        while (iterationIntratype < intratypeIterations) {
            #pragma omp parallel for
            for(uint i = 0; i < typesFiltered.size(); i++){
                std::vector<std::string> nodeNames = typeToNodeNames[i];
                std::cout << "[LOG] computation of perturbation for iteration intertype-intratype ("+ std::to_string(iterationIntertype) + "<->"+ std::to_string(iterationIntratype) + ") for type (" + types[i]<<std::endl; 
                
                if (saturation) {
                    if(vm.count("saturationTerm") == 0){
                        std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced2((iterationIntertype*intratypeIterations + iterationIntratype)*timestep, saturation = true); // TODO check if iteration intratype should be multiplied by iteration intertype
                    } else if (vm.count("saturationTerm") >= 1) {
                        //TODO create saturation vector
                        double saturationTerm = vm["saturationTerm"].as<double>();
                        //TODO TEST
                        std::vector<double> saturationVector = std::vector<double>(graphsNodes[invertedTypesIndexes[i]].size(),saturationTerm);
                        std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced2((iterationIntertype*intratypeIterations + iterationIntratype)*timestep, saturation = true, saturationVector); // TODO check if iteration intratype should be multiplied by iteration intertype
                    }
                } else{
                    std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced2(iterationIntertype*intratypeIterations + iterationIntratype, saturation = false); // TODO check if iteration intratype should be multiplied by iteration intertype
                }
            }
            //save output values
            for(uint i = 0; i < typesFiltered.size(); i++){
                std::vector<std::string> nodeNames = typeToNodeNames[i];
                //TODO change how to save files to get more information about intratype and intertype iterations
                saveNodeValues(outputFoldername, iterationIntertype*intratypeIterations + iterationIntratype, types[i], typeComputations[i]->getOutputAugmented(), nodeNames,ensembleGeneNames);
            }
            // std::cout<< "[DEBUG] output values before updating input"<<std::endl;
            // for(uint i = 0; i < types.size(); i++){
            //     std::cout << "[DEBUG] type " << types[i] << " values: ";
            //     for(uint j = 0; j < typeComputations[i]->getOutputAugmented().size(); j++){
            //         std::cout << typeComputations[i]->getOutputAugmented()[j] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            //update input
            for(uint i = 0; i < typesFiltered.size(); i++){
                //If conservation of the initial values is required, the input is first updated with the initial norm value
                if (conservateInitialNorm) {
                    std::vector<double> inputInitial = inputInitials[i];
                    double initialNorm = vectorNorm(inputInitial);
                    double outputNorm = vectorNorm(typeComputations[i]->getOutputAugmented());
                    double normRatio = initialNorm/outputNorm;
                    std::vector<double> newInput = vectorScalarMultiplication(typeComputations[i]->getOutputAugmented(),normRatio);
                    std::cout << "[LOG] update input with conservation of the initial perturbation for iteration intertype-intratype ("+ std::to_string(iterationIntertype) + "<->"+ std::to_string(iterationIntratype) + ") for type (" + types[i]<<std::endl;
                    typeComputations[i]->updateInput(newInput,true);
                } else {
                    std::cout << "[LOG] update input for iteration intertype-intratype ("+ std::to_string(iterationIntertype) + "<->"+ std::to_string(iterationIntratype) + ") for type (" + types[i]<<std::endl;
                    typeComputations[i]->updateInput(std::vector<double>(),true);
                }
                
            }
            // std::cout<< "[DEBUG] input values after updating input"<<std::endl;
            // for(uint i = 0; i < types.size(); i++){
            //     std::cout << "[DEBUG] type " << types[i] << " values: ";
            //     for(uint j = 0; j < typeComputations[i]->getInputAugmented().size(); j++){
            //         std::cout << typeComputations[i]->getInputAugmented()[j] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            iterationIntratype++;
        }
        //update input with virtual node values update
        
        // std::cout<< "[DEBUG] input values before updating with virtual"<<std::endl;
        // for(uint i = 0; i < types.size(); i++){
        //     std::cout << "[DEBUG] type " << types[i] << " values: ";
        //     for(uint j = 0; j < typeComputations[i]->getInputAugmented().size(); j++){
        //         std::cout << typeComputations[i]->getInputAugmented()[j] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        for (uint i = 0; i < typesFiltered.size(); i++) {
            //queuesTypeTypes[i] = typeComputations[i]->computeAugmentedPerturbation();
            //TODO when computation will be done in parallel, this step should wait for all the computations of the other adjacent types to finish
            //also take into account REpast framework
            for(uint j = 0; j < typesFiltered.size(); j++){
                if(i==j){
                    if(sameTypeCommunication) typeComputations[i]->setInputVinForCell(typesFiltered[j], typeComputations[j]->getVirtualOutputForCell(typesFiltered[i]));
                } else {
                    typeComputations[i]->setInputVinForCell(typesFiltered[j], typeComputations[j]->getVirtualOutputForCell(typesFiltered[i]));
                }
            }
        }
        // std::cout<< "[DEBUG] input values after updating with virtual"<<std::endl;
        // for(uint i = 0; i < types.size(); i++){
        //     std::cout << "[DEBUG] type " << types[i] << " values: ";
        //     for(uint j = 0; j < typeComputations[i]->getInputAugmented().size(); j++){
        //         std::cout << typeComputations[i]->getInputAugmented()[j] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        
        iterationIntertype++;

    }
    


    //cleaning memory
    // for(uint i = 0; i< types.size(); i++){
    //     delete typeComputations[i];
    // }
    // delete [] typeComputations;
    //delete graph
    

    
    return 0;
}