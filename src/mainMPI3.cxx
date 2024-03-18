#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <boost/program_options/value_semantic.hpp>
#include <iostream>
#include <boost/program_options.hpp>
#include <map>
#include <sys/types.h>
#include <tuple>
#include "Computation.h"
#include "PropagationModel.hxx"
#include "PropagationModelOriginal.hxx"
#include "PropagationModelNeighbors.hxx"
#include "PropagationModelCustom.hxx"
#include "ConservationModel.h"
#include "DissipationModel.h"
#include "DissipationModelPow.h"
#include "DissipationModelRandom.h"
#include "DissipationModelScaled.h"
#include "WeightedEdgeGraph.h"
#include "utilities.h"
#include "CustomScalingFunctions.h"
#include "Logger.hxx"



int main(int argc, char** argv) {    
    //program options
    bool ensembleGeneNames=false;
    bool sameTypeCommunication=false;
    bool saturation=false;
    bool conservateInitialNorm=false;
    bool undirected = false;
    bool undirectedTypeEdges = false;
    bool resetVirtualOutputs = false;
    std::string logMode="";
    std::string virtualNodesGranularity = "type";
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    std::string performanceFilename = "";
    desc.add_options()
        ("help", "() print help section")//<initialPerturbationPerType>.tsv [<subtypes>.txt] [<typesInteraction>.tsv]\nFILE STRUCTURE SCHEMA:\ngraph.tsv\nstart end weight\n<gene1> <gene2>  <0.something>\n...\n\n\ninitialPerturbationPerType.tsv\n type1 type2 ... typeN\ngene1 <lfc_type1:gene1> <lfc_type2:gene1> ... <lfc_typeN:gene1>\ngene1 <lfc_type1:gene2> <lfc_type2:gene2> ... <lfc_typeN:gene2>\n...\n\n\ntypesInteraction.tsv\nstartType:geneLigand endType:geneReceptor weight\n<type1:geneLigand> <type2:genereceptor>  <0.something>\n...\n\n\nsubtypes.txt\ntype1\ntype3\n...")
        ("fUniqueGraph", po::value<std::string>(), "(string) graph filename, for an example graph see in resources. NOTE: if this option is chosen, graphsFilesFolder cannot be used. For an example see in data data/testdata/testGraph/edges-Graph1-general.tsv")
        ("fInitialPerturbationPerType", po::value<std::string>(), "(string) initialPerturbationPerType matrix filename, for an example see in data data/testdata/testGraph/initialValues-general.tsv")
        ("subtypes", po::value<std::string>(), "subtypes filename, for an example see in data, see data/testdata/testGraph/subcelltypes.txt")
        ("initialPerturbationPerTypeFolder", po::value<std::string>(), "(string) initialPerturbationPerType folder, for an example see in data data/testdata/testGraph/initialValues")
        ("typeInteractionFolder", po::value<std::string>(), "(string) directory for the type interactions, for an example see in data data/testdata/testHeterogeneousGraph/interactions")
        ("nodeDescriptionFile", po::value<std::string>(), "(string) node description file, used to generate the output description, if not specified no names are used. for an example see in data resources/graphs/metapathwayNew/nodes.tsv")
        ("ensembleGeneNames",po::bool_switch(&ensembleGeneNames),"() use ensemble gene names, since the graph used in resources have entrez_ids, a map will be done from ensemble to entrez, the map is available in resources")
        ("sameTypeCommunication",po::bool_switch(&sameTypeCommunication),"() use same type communication, since it is not permitted as the standard definition of the model, this adds a virtual node for the same type type")
        ("outputFolder",po::value<std::string>()->required(),"(string) output folder for output of the algorithm at each iteration")
        ("intertypeIterations",po::value<uint>(),"(positive integer) number of iterations for intertype communication")
        ("intratypeIterations",po::value<uint>(),"(positive integer) number of iterations for intratype communication")
        ("timestep",po::value<double>(),"timestep to use for the iteration, the final time is iterationIntracell*iterationIntercell*timestep")
        ("dissipationModel",po::value<std::string>(),"(string) the dissipation model for the computation, available models are: 'none (default)','power','random','periodic','scaled' and 'custom'")
        ("dissipationModelParameters",po::value<std::vector<double>>()->multitoken(),"(string) the parameters for the dissipation model, for the power dissipation indicate the base, for the random dissipation indicate the min and max value, for the periodic dissipation indicate the period")
        ("graphsFilesFolder",po::value<std::string>(),"(string) graphs (pathways or other types of graphs) file folder, for an example see in data data/testdata/testHeterogeneousGraph/graphsDifferentStructure")
        ("conservationModel",po::value<std::string>(),"(string) the conservation model used for the computation, available models are: 'none (default)','scaled','random' and 'custom' ")
        ("conservationModelParameters", po::value<std::vector<double>>()->multitoken(),"(vector<double>) the parameters for the dissipation model, for the scaled parameter the constant used to scale the conservation final results, in the case of random the upper and lower limit (between 0 and 1)")
        ("propagationModel",po::value<std::string>(),"(string) the propagation model used for the computation, available models are: 'default(pseudoinverse creation)','scaled (pseudoinverse * scale parameter)', neighbors(propagate the values only on neighbors at every iteration and scale parameter) and 'customScaling' (pseudoinverse*scalingFunction(parameters)), 'customScalingNeighbors' (neighbors propagation and scalingFunction(parameters)), 'customPropagation' (custom scaling function and custom propagation function defined in src/PropagationModelCustom) ")
        ("propagationModelParameters", po::value<std::vector<double>>()->multitoken(),"(vector<double>) the parameters for the propagation model, for the scaled parameter the constant used to scale the conservation final results")
        ("saturation",po::bool_switch(&saturation),"use saturation of values, default to 1, if another value is needed, use the saturationTerm")
        ("saturationTerm",po::value<double>(),"defines the limits of the saturation [-saturationTerm,saturationTerm]")
        ("conservateInitialNorm",po::bool_switch(&conservateInitialNorm), "conservate the initial euclidean norm of the perturbation values, that is ||Pn|| <= ||Initial||, default to false")
        ("undirectedEdges",po::bool_switch(&undirected), "edges in the graphs are undirected")
        ("undirectedTypeEdges",po::bool_switch(&undirectedTypeEdges), "edges between types are undirected")
        ("resetVirtualOutputs",po::bool_switch(&resetVirtualOutputs), "reset the virtual outputs to 0 at each iteration, default to false")
        ("virtualNodesGranularity", po::value<std::string>(), "(string) granularity of the virtual nodes, available options are: 'type', 'node'(unstable), 'typeAndNode', default to type")
        ("virtualNodesGranularityParameters", po::value<std::vector<std::string>>()->multitoken(), "(vector<string>) parameters for the virtual nodes granularity, NOT USED for now")
        ("loggingOptions",po::value<std::string>(&logMode),"(string) logging options, available options are: 'all','none', default to all")
        ("savePerformance",po::value<std::string>(&performanceFilename), "(string) output performance (running time, number of total nodes, number of communities, number of total edges) to the defined file, if nothing is specified the performance are not saved")
    ;
    //TODO add additional boolean parameter to control if the graph names are not genes and the algorithm should use the graph names directly, no conversion or mapping

    

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    std::string filename,subtypesFilename,typesFilename,typesInteractionFoldername,typesInitialPerturbationMatrixFilename,graphsFilesFolder, typeInitialPerturbationFolderFilename,outputFoldername;
    int intertypeIterations,intratypeIterations;
    DissipationModel* dissipationModel = nullptr;
    ConservationModel* conservationModel = nullptr;
    double timestep = 1;


    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    //logging options
    Logger logger(std::cout);
    if(vm.count("loggingOptions")){
        if(logMode == "all"){
            std::cout << "[LOG] logging options set to all"<<std::endl;
            logger.enable();
        } else if (logMode == "none"){
            std::cout << "[LOG] logging options set to none"<<std::endl;
            logger.disable();
        } else {
            std::cout << "[LOG] logging options set to default (all)"<<std::endl;
            logger.enable();
        }
    } else {
        std::cout << "[LOG] logging options set to default (all)"<<std::endl;
        logger.enable();
    }


    //controls over impossible configurations
    if(vm.count("fUniqueGraph") == 0 && vm.count("graphsFilesFolder") == 0){
        //no unique graph of folder of the graphs was set
        std::cerr << "[ERROR] no unique graph filename or folder was set to get the graphs, set one "<<std::endl;
        return 1;
    }

    if(vm.count("fInitialPerturbationPerType") == 0 && vm.count("initialPerturbationPerTypeFolder") == 0){
        //no way of getting the initial perturbation values
        std::cerr << "[ERROR] no matrix for the initial values was passed as filename or single vector in files contained in the folder specified was set, set one "<<std::endl;
        return 1;
    }

    if(vm.count("fInitialPerturbationPerType") && vm.count("graphsFilesFolder")){
        //unstable configuration of different graphs and single matrix with the same nodes
        logger << "[WARNING] unstable configuration of different graphs and a single matrix with the initial perturbations"<<std::endl;
    }

    if(saturation && conservateInitialNorm){
        std::cerr << "[ERROR] saturation and conservateInitialNorm cannot be both true, aborting"<<std::endl;
        return 1;
    }


    if(vm.count("graphsFilesFolder") && vm.count("fUniqueGraph")){
        std::cerr << "[ERROR] fUniqueGraph and graphsFilesFolder were both set. Aborting\n";
        return 1;
    }
    if(vm.count("initialPerturbationPerTypeFolder") && vm.count("fInitialPerturbationPerType")){
        std::cerr << "[ERROR] fInitialPerturbationPerType and initialPerturbationPerTypeFolder were both set. Aborting\n";
        return 1;
    }

    // reading the parameters

    if (vm.count("intertypeIterations")) {
        logger << "[LOG] iterations intertype set to " 
    << vm["intertypeIterations"].as<uint>() << ".\n";
        intertypeIterations = vm["intertypeIterations"].as<uint>();
    } else {
        logger << "[LOG] iterations intertype not set, set to default: 10 iterations \n";
        intertypeIterations = 10;
    }

    if (vm.count("intratypeIterations")) {
        logger << "[LOG] iterations intratype set to " 
    << vm["intratypeIterations"].as<uint>() << ".\n";
        intratypeIterations = vm["intratypeIterations"].as<uint>();
    } else {
        logger << "[LOG] iterations intratype not set, set to default: 5 iterations \n";
        intratypeIterations = 5;
    }

    // reading granularity
    if(vm.count("virtualNodesGranularity")){
        virtualNodesGranularity = vm["virtualNodesGranularity"].as<std::string>();
        logger << "[LOG] virtual nodes granularity set to " << virtualNodesGranularity << std::endl;
        // controls over the value
        if(virtualNodesGranularity != "type" && virtualNodesGranularity != "node" && virtualNodesGranularity != "typeAndNode"){
            std::cerr << "[ERROR] virtual nodes granularity must be one of the following: 'type', 'node' or 'typeAndNode': aborting"<<std::endl;
            return 1;
        }
        if(virtualNodesGranularity == "node"){
            std::cerr << "[WARNING] virtual nodes granularity set to 'node', this option is unstable and not fully implemented: aborting"<<std::endl;
            return 1;
        }
    } else {
        logger << "[LOG] virtual nodes granularity not set, set to default: type \n";
        virtualNodesGranularity = "type";
    }

    // reading granularity parameters (not used for now)
    if(vm.count("virtualNodesGranularityParameters")){
        std::vector<std::string> virtualNodesGranularityParameters = vm["virtualNodesGranularityParameters"].as<std::vector<std::string>>();
        logger << "[LOG] virtual nodes granularity parameters set "<< std::endl;
        logger << "[WARNING] virtual nodes granularity parameters are not used for now"<<std::endl;
    } else {
        // logger << "[LOG] virtual nodes granularity parameters not set, set to default: empty vector \n";  // not used for now
    }

    //logging timestep settings
    if(vm.count("timestep")){
        logger << "[LOG] timestep set to " 
    << vm["timestep"].as<double>() << ".\n";
        timestep = vm["timestep"].as<double>();
    } else {
        logger << "[LOG] timestep not set, set to default (1)"<<std::endl;
    }

    //logging edges direction in the graphs
    if(undirected){
        logger << "[LOG] undirectedEdges specified, undirected edges in the graphs"<<std::endl;
    } else {
        logger << "[LOG] undirectedEdges not specified, directed edges in the graphs(only the edges specified in the graph files will be added)"<<std::endl;
    }

    //logging edges direction in the graphs
    if(undirectedTypeEdges){
        logger << "[LOG] undirectedTypeEdges specified, undirected edges between types"<<std::endl;
    } else {
        logger << "[LOG] undirectedTypeEdges not specified, directed edges between types"<<std::endl;
    }

    //logging virtual nodes reset

    if(resetVirtualOutputs){
        logger << "[LOG] resetVirtualOutputs specified, virtual outputs will be reset to 0 after each inter-propagation"<<std::endl;
    } else {
        logger << "[LOG] resetVirtualOutputs not specified, virtual outputs will not be reset to 0 after each inter-propagation"<<std::endl;
    }


    if (vm.count("fUniqueGraph")) {
        logger << "[LOG] file for the graph was set to " 
    << vm["fUniqueGraph"].as<std::string>() << ".\n";
        filename = vm["fUniqueGraph"].as<std::string>();
        if(!fileExistsPath(filename)){
            std::cerr << "[ERROR] file for the graph do not exist: aborting"<<std::endl;
            return 1;
        }
    } else if(vm.count("graphsFilesFolder")){
        logger << "[LOG] folder for the graphs was set to " 
    << vm["graphsFilesFolder"].as<std::string>() << ".\n";
        graphsFilesFolder = vm["graphsFilesFolder"].as<std::string>();
        if(!folderExists(graphsFilesFolder)){
            std::cerr << "[ERROR] folder for the graphs do not exist: aborting"<<std::endl;
            return 1;
        }
    }
    if (vm.count("fInitialPerturbationPerType")) {
        logger << "[LOG] file for the initialPerturbationPerType matrix was set to " 
    << vm["fInitialPerturbationPerType"].as<std::string>() << ".\n";
        typesInitialPerturbationMatrixFilename = vm["fInitialPerturbationPerType"].as<std::string>();
        if(!fileExistsPath(typesInitialPerturbationMatrixFilename)){
            std::cerr << "[ERROR] file for the initialPerturbationPerType does not exist: aborting"<<std::endl;
            return 1;
        }
    } else if (vm.count("initialPerturbationPerTypeFolder")) {
        logger << "[LOG] folder for the initialPerturbationPerType was set to "
    << vm["initialPerturbationPerTypeFolder"].as<std::string>() << ".\n";
        typeInitialPerturbationFolderFilename = vm["initialPerturbationPerTypeFolder"].as<std::string>();
        if(!folderExists(typeInitialPerturbationFolderFilename)){
            std::cerr << "[ERROR] folder for the initialPerturbationPerType do not exist: aborting"<<std::endl;
            return 1;
        }
    }

    if (vm.count("typeInteractionFolder")) {
        logger << "[LOG] folder for the type interactions was set to " 
    << vm["typeInteractionFolder"].as<std::string>() << ".\n";
        typesInteractionFoldername = vm["typeInteractionFolder"].as<std::string>();
        if(!folderExists(typesInteractionFoldername)){
            std::cerr << "[ERROR] folder for the type interactions do not exist: aborting"<<std::endl;
            return 1;
        }
    } else {
        logger << "[LOG] typeInteractionFolder folder was not set. computing without taking into account type interactions\n";
    }
    if (vm.count("outputFolder")) {
        logger << "[LOG] output folder  was set to " 
    << vm["outputFolder"].as<std::string>() << ".\n";
        outputFoldername = vm["outputFolder"].as<std::string>();
        if(!folderExists(outputFoldername)){
            std::cerr << "[WARNING] folder for the output do not exist: creating the folder"<<std::endl;
            if(!createFolder(outputFoldername)){
                std::cerr << "[ERROR] folder for the output could not be created: aborting"<<std::endl;
                return 1;
            }
        }
    } else {
        std::cerr << "[ERROR] output folder was not set. aborting\n";
        return 1;
    }
    if (vm.count("dissipationModel")) {
        logger << "[LOG] dissipation model was set to "
    << vm["dissipationModel"].as<std::string>() << ".\n";
        std::string dissipationModelName = vm["dissipationModel"].as<std::string>();
        if(dissipationModelName == "none"){
            logger << "[LOG] dissipation model set to default (none)\n";
            dissipationModel = new DissipationModelScaled([](double time)->double{return 0;});
        } else if(dissipationModelName == "power"){
            if (vm.count("dissipationModelParameters")) {
                logger << "[LOG] dissipation model parameters for power dissipation were declared to be" << vm["dissipationModelParameters"].as<std::vector<double>>()[0] << ".\n";
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
                logger << "[LOG] dissipation model parameters were declared to be "
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
                logger << "[LOG] dissipation model parameters were declared to be "
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
                    logger << "[LOG] dissipation model parameters were set to Amplitude:"
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
            logger << "[LOG] dissipation model was set to custom, if the function is not correctly defined there will be errors\n " << std::endl;
            dissipationModel = new DissipationModelScaled(getDissipationScalingFunction());
        } else {
            std::cerr << "[ERROR] dissipation model scale function is not any of the types. Conservation model scale functions available are none(default), scaled, random and custom \n";
            return 1;
        }
    } else { //dissipation model set to default (none)
        logger << "[LOG] dissipation model was not set. set to default (none)\n";
        dissipationModel = new DissipationModelScaled([](double time)->double{return 0;});
    }


    if (vm.count("conservationModel")) {
        logger << "[LOG] conservation model was set to "
    << vm["conservationModel"].as<std::string>() << ".\n";
        std::string conservationModelName = vm["conservationModel"].as<std::string>();
        if(conservationModelName == "none"){
            logger << "[LOG] conservation model set to default (none)\n";
            conservationModel = new ConservationModel([](double time)->double{return 0;});
        } else if (conservationModelName == "scaled"){
            if (vm.count("conservationModelParameters")) {
                logger << "[LOG] conservation model parameters were declared to be "
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
                logger << "[LOG] conservation model parameters were declared to be "
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
            logger << "[LOG] conservation model was set to custom, if the custom function defined for scaling is not correctly implemented, there will be errors\n " << std::endl;
            conservationModel = new ConservationModel(getConservationScalingFunction());
        } else {
            std::cerr << "[ERROR] conservation model scale function is not any of the types. Conservation model scale functions available are none(default), scaled, random and custom \n";
            return 1;
        }
    } else {
        logger << "[LOG] conservation model was not set. set to default (none)\n";
        conservationModel = new ConservationModel([](double time)->double{return 0;});
    }

    //logging if saturation is set and saturation parameters are set
    if (saturation) {
        if(vm.count("saturationTerm") == 0){
            logger << "[LOG] saturation term not specified, using the interval [-1,1]"<<std::endl;
        } else if(vm.count("saturationTerm") == 1){
            double saturationTerm = vm["saturationTerm"].as<double>();
            logger << "[LOG] saturation term specified, using the interval [-" << saturationTerm << "," << saturationTerm << "]"<<std::endl;
        } else {
            std::cerr << "[ERROR] saturation term specified more than once, possibility of using more values not yet implemented: aborting"<<std::endl;
            return 1;
        }
    }
    //end program options section

    std::string nodesDescriptionFilename="";
    if(ensembleGeneNames){
        logger <<"[LOG] mapping ensemble gene names to entrez ids"<<std::endl;
    } else if(vm.count("nodeDescriptionFile")){
        logger <<"[LOG] using node description file to get the names of the nodes in the graphs"<<std::endl;
        nodesDescriptionFilename = vm["nodeDescriptionFile"].as<std::string>();
    } else {
        logger <<"[LOG] no nodes description"<<std::endl;
    }


    // these types will be used for indexing the single processes, rank 0 is the master process, rank 1 will get the first type, rank 2 the second and so on
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
        logger << "[LOG] subtypes filename set to "
    << vm["subtypes"].as<std::string>() << ".\n";
        subtypesFilename = vm["subtypes"].as<std::string>();
        subtypes = getVectorFromFile<std::string>(subtypesFilename);
    }else{
        logger << "[LOG] subtypes filename not set, set to default: all types \n";
        subtypes = types;
    }
    
    //filter types with the subtypesresetVirtualOutputs
    types = vectorsIntersection(types, subtypes);
    if (types.size() == 0) {
        std::cerr << "[ERROR] no types in common between the types and subtypes: aborting"<<std::endl;
        return 1;
    }


    // take starting time before initializing MPI and the computation
    auto start = std::chrono::steady_clock::now();
    // initialize MPI
    MPI_Init(&argc, &argv);

    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // control if the number of processes is greater than the number of types, exit if true (useless process are not accepted)
    if (numProcesses > SizeToInt(types.size())) {
        std::cerr << "[ERROR] number of processes is greater than the number of types: aborting"<<std::endl;
        std::cerr << "[ERROR] number of processes: " << numProcesses << std::endl;
        std::cerr << "[ERROR] number of types: " << types.size() << std::endl;
        return 1;
    }

    // workloads for each process, 
    int workloadPerProcess = round((types.size() + 0.0) / numProcesses); // one process minimum for each type, if number of processes = number of types then one process for each type, rounded to account for close number of types and processes(almost 1 process for every type) and to not leave a process without work
    int startIdx = rank * workloadPerProcess;
    int endIdx = (rank == numProcesses - 1) ? types.size() : (rank + 1) * workloadPerProcess;
    
    int finalWorkload = endIdx - startIdx;

    //map types to rank
    std::map<std::string, int> typeToRank;
    for (int i = 0; i < SizeToInt(types.size()); ++i) {
        int rankType;
        if (i >= (numProcesses - 1)*workloadPerProcess) {
            rankType = numProcesses - 1;
        } else {
            rankType = i / workloadPerProcess;
        }
        typeToRank[types[i]] = rankType;
    }

    //use the number of types for workload to allocate an array of pointers to contain the graph for each type
    WeightedEdgeGraph **graphs = new WeightedEdgeGraph*[finalWorkload];
    std::vector<std::vector<std::string>> graphsNodes;
    std::vector<std::vector<std::string>> graphsNodesAll; // used only initially to read the values, contains all types
    std::unordered_map<std::string, std::vector<std::string>> typeToNodeNamesMap; // map from all types to the node names, not only the ones in the workload, no virtual nodes
    std::vector<std::pair<std::vector<std::string>,std::vector<std::tuple<std::string,std::string,double>>>> namesAndEdges;
    // a single graph is used for all the types
    if(vm.count("fUniqueGraph")){
        namesAndEdges.push_back(edgesFileToEdgesListAndNodesByName(filename));
        graphsNodes.push_back(namesAndEdges[0].first);
        graphs[0] = new WeightedEdgeGraph(graphsNodes[0]);
        for(int i = 1; i < finalWorkload; i++){
            namesAndEdges.push_back(namesAndEdges[0]);
            graphsNodes.push_back(namesAndEdges[0].first);
            graphs[i] = graphs[0];
        }
        // create the map also for the use case of a single graph, TODO only use the map for the single graph case without filling all the maps
        for(uint i = 0; i < types.size(); i++){
            typeToNodeNamesMap[types[i]] = namesAndEdges[0].first;
        }

    } else if (vm.count("graphsFilesFolder")) { // the graphs are in a folder, each graph is a type
        auto allGraphs = edgesFileToEdgesListAndNodesByNameFromFolder(graphsFilesFolder);
        auto typesFromFolder = allGraphs.first;
        if(typesFromFolder.size() != types.size()){
            std::cerr << "[ERROR] types from folder and types from file do not match: aborting"<<std::endl;
            return 1;
        }
        for (uint i = 0; i<typesFromFolder.size(); i++){ //TODO map the types from the folder to the types from the file
            if(typesFromFolder[i] != types[i]){
                std::cerr << "[ERROR] types from folder and types from file do not match: aborting"<<std::endl;
                return 1;
            }
        }
        namesAndEdges = allGraphs.second;
        for(uint i = 0; i < types.size(); i++){
            graphsNodesAll.push_back(namesAndEdges[i].first);
            typeToNodeNamesMap[types[i]] = namesAndEdges[i].first;
        }
        for(int i = startIdx; i < endIdx; i++){
            graphsNodes.push_back(namesAndEdges[i].first);
            graphs[i-startIdx] = new WeightedEdgeGraph(graphsNodes[i-startIdx]);
        }
    } 

    //add the edges to the graphs
    if(vm.count("fUniqueGraph")){
        for(auto edge = namesAndEdges[0].second.cbegin() ; edge != namesAndEdges[0].second.cend(); edge++ ){
            graphs[0]->addEdge(std::get<0> (*edge), std::get<1> (*edge) ,std::get<2>(*edge) ,!undirected);
        }
    } else if (vm.count("graphsFilesFolder")) {
        for(int i = startIdx; i < endIdx; i++){
            for(auto edge = namesAndEdges[i].second.cbegin() ; edge != namesAndEdges[i].second.cend(); edge++ ){
                graphs[i-startIdx]->addEdge(std::get<0> (*edge), std::get<1> (*edge) ,std::get<2>(*edge) ,!undirected);
            }
        }
    }

    //get initial input values
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::vector<double>>> initialValues;
    std::vector<std::vector<double>> inputInitials;
    if(vm.count("fInitialPerturbationPerType")){
        logger << "[LOG] initial perturbation per type specified, using the file "<<typesInitialPerturbationMatrixFilename<<std::endl;
        initialValues = logFoldChangeMatrixToCellVectors(typesInitialPerturbationMatrixFilename,graphsNodes[0],subtypes,ensembleGeneNames);
    } else if (vm.count("initialPerturbationPerTypeFolder")){
        logger << "[LOG] initial perturbation per type specified, using the folder "<<typeInitialPerturbationFolderFilename<<std::endl;
        initialValues = logFoldChangeCellVectorsFromFolder(typeInitialPerturbationFolderFilename,types,graphsNodesAll,subtypes,ensembleGeneNames);
    } else {
        std::cerr << "[ERROR] no initial perturbation file or folder specified: aborting"<<std::endl;
        return 1;
    }
    std::vector<std::string> initialNames = std::get<0>(initialValues);
    inputInitials = std::get<2>(initialValues);
    std::vector<std::string> typesFromValues = std::get<1>(initialValues);
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

    // TODO get rid of the augment graph function and use the addition of nodes or edges directly
    Computation** typeComputations = new Computation*[finalWorkload];
    int indexComputation = 0;
    std::vector<int> typesIndexes = std::vector<int>(finalWorkload,-1); 
    std::vector<int> invertedTypesIndexes = std::vector<int>(finalWorkload,-1); 
    for(int i = 0; i < finalWorkload; i++){
        if(indexMapGraphTypesToValuesTypes[i+startIdx] == -1){
            logger << "[LOG] type "<<types[i+startIdx]<<" not found in the initial perturbation files, using zero vector as input"<<std::endl;
            std::vector<double> input = std::vector<double>(graphsNodes[i].size(),0);
            Computation* tmpCompPointer = new Computation(types[i+startIdx],input,graphs[i],graphsNodes[i]);   
            tmpCompPointer->setDissipationModel(dissipationModel);
            tmpCompPointer->setConservationModel(conservationModel);
            typeComputations[indexComputation] = tmpCompPointer;
            //No inverse computation with the augmented graph since virtual nodes edges are not yet inserted
            // TODO generalize by removing the type granularity in this code, that is by considering only the types that are encoded?
            if(virtualNodesGranularity == "type"){
                typeComputations[indexComputation]->augmentGraphNoComputeInverse(types,std::vector<std::pair<std::string,std::string>>(),std::vector<double>(), true); //self included since the code in MPI needs it
            } else if (virtualNodesGranularity == "typeAndNode"){
                typeComputations[indexComputation]->augmentGraphNoComputeInverse(std::vector<std::string>(), std::vector<std::pair<std::string,std::string>>(), std::vector<double>(), false); //no types are passed since the virtual nodes will be added to the graph in the interaction section of this code
            }
        } else {
            int index = indexMapGraphTypesToValuesTypes[i+startIdx];
            std::vector<double> input = inputInitials[index];
            Computation* tmpCompPointer = new Computation(types[i+startIdx],input,graphs[i],graphsNodes[i]); 
            tmpCompPointer->setDissipationModel(dissipationModel);
            tmpCompPointer->setConservationModel(conservationModel);
            typeComputations[indexComputation] = tmpCompPointer;
            //No inverse computation with the augmented graph since virtual nodes edges are not yet inserted
            if(virtualNodesGranularity == "type"){
                typeComputations[indexComputation]->augmentGraphNoComputeInverse(types,std::vector<std::pair<std::string,std::string>>(),std::vector<double>(), true); //self included since the code in MPI needs it
            } else if (virtualNodesGranularity == "typeAndNode"){
                typeComputations[indexComputation]->augmentGraphNoComputeInverse(std::vector<std::string>(), std::vector<std::pair<std::string,std::string>>(), std::vector<double>(), false); //no types are passed since the virtual nodes will be added to the graph in the interaction section of this code
            }
        }
        typesIndexes[i] = indexComputation;
        invertedTypesIndexes[indexComputation] = i;
        indexComputation++;

    }

    auto allFilesInteraction = get_all(typesInteractionFoldername,".tsv");
    // define the map for the type interactions, an hash function should be defined for the pair of strings used as the identifier of the interaction
    std::unordered_map<std::pair<std::string, std::string>, std::unordered_set<int>, hash_pair_strings> interactionBetweenTypesMap;
    std::unordered_map<std::tuple<std::string, std::string, std::string, std::string>, std::unordered_set<int>, hash_quadruple_strings> interactionBetweenTypesFinerMap;
    for(auto typeInteractionFilename = allFilesInteraction.cbegin() ; typeInteractionFilename != allFilesInteraction.cend() ; typeInteractionFilename++){
        std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::unordered_set<int>,double>>> typeInteractionsEdges;
        if (subtypes.size() == 0) {
            // TODO add different contact times inside the network (quite difficult since the structure of the graphs is static)
            // the above can be implemented with the use of different matrices for every single type(from 1 to max(contactTimes)), where the matrix represent the current state of the network
            // another possibility is to use two matrices for every type-agent, one for the whole network without contact times and one is used to store the current network state at iteration i 
            // SOLUTION: granularity
            typeInteractionsEdges  = interactionContactsFileToEdgesListAndNodesByName(*typeInteractionFilename, types, intertypeIterations, ensembleGeneNames, virtualNodesGranularity, typeToNodeNamesMap, undirectedTypeEdges);
        } else {
            typeInteractionsEdges = interactionContactsFileToEdgesListAndNodesByName(*typeInteractionFilename, subtypes, intertypeIterations, ensembleGeneNames, virtualNodesGranularity, typeToNodeNamesMap, undirectedTypeEdges);
        }
        #pragma omp parallel for
        for (int i = 0; i < finalWorkload;i++) {
            if(typeInteractionsEdges.first.contains(types[i+startIdx]) && typesIndexes[i] != -1){
                // granularity is already considered in the function that reads from the file previously called
                typeComputations[typesIndexes[i]]->addEdgesAndNodes(typeInteractionsEdges.first[types[i+startIdx]], false, false); // no inverse computation since it is done in the propagation model
            }
        }
        for(auto edge = typeInteractionsEdges.second.cbegin() ; edge != typeInteractionsEdges.second.cend(); edge++ ){
            // first two types are the nodes in the two networks/types ,types are the third and fourth element of the tuple, while the fifth is the set of contact times
            // startNodeName, endNodeName, startType, endType, contactTimes, weight
            std::pair<std::string,std::string> keyTypes = std::make_pair(std::get<2> (*edge), std::get<3> (*edge));
            std::tuple<std::string,std::string,std::string,std::string> keyTypesFiner = std::make_tuple(std::get<0> (*edge), std::get<1> (*edge), std::get<2> (*edge), std::get<3> (*edge));
            
            // TESTING
            // printing keytypesFiner
            // std::cout << "keyTypesFiner: (" << std::get<0> (keyTypesFiner) << " " << std::get<1> (keyTypesFiner) << " " << std::get<2> (keyTypesFiner) << " " << std::get<3> (keyTypesFiner)<< ")" << std::endl;
            // TESTING

            if(interactionBetweenTypesMap.contains(keyTypes)){
                interactionBetweenTypesMap[keyTypes].insert(std::get<4>(*edge).begin(),std::get<4>(*edge).end()); // directly inserting means the union of the two sets
            } else {
                interactionBetweenTypesMap[keyTypes] = std::get<4>(*edge);
            }
            
            if(interactionBetweenTypesFinerMap.contains(keyTypesFiner)){
                interactionBetweenTypesFinerMap[keyTypesFiner].insert(std::get<4>(*edge).begin(),std::get<4>(*edge).end()); // directly inserting means the union of the two sets
            } else {
                interactionBetweenTypesFinerMap[keyTypesFiner] = std::get<4>(*edge);
            }
            
        }
    }

    // create a map that maps couples of strings (source type and target type) to a vector of pairs of strings, representing how the virtual outputs are mapped in the subarray passed to MPI send 
    std::unordered_map<std::pair<std::string, std::string>, std::vector<std::pair<std::string, std::string>>,hash_pair_strings> typesPairMappedVirtualOutputsVectors;
    // create a map that maps couples of strings (source type and target type) to a vector of pairs of strings, representing how the virtual inputs are mapped in the subarray passed to MPI send
    std::unordered_map<std::pair<std::string, std::string>, std::vector<std::pair<std::string, std::string>>,hash_pair_strings> typesPairMappedVirtualInputsVectors;

    // create a map that maps couples of integers (source rank and target rank) to a vector of pairs of strings(virtual output name, virtual input name), representing how the virtual nodes are mapped in the array passed to MPI send
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<std::string, std::string>>,hash_pair_ints> ranksPairMappedVirtualNodesVectors;

    
    // populate the maps
    for(auto interaction = interactionBetweenTypesFinerMap.cbegin() ; interaction != interactionBetweenTypesFinerMap.cend(); interaction++ ){
        std::string startNodeName = std::get<0> (interaction->first);
        std::string endNodeName = std::get<1> (interaction->first);
        std::string startType = std::get<2> (interaction->first);
        std::string endType = std::get<3> (interaction->first);
        std::pair<std::string, std::string> keyTypes = std::make_pair(startType, endType);
        // inverted interaction only used when specified by the command line argument (undirectedTypeEdges)
        std::pair<std::string, std::string> keyTypesInverted = std::make_pair(endType, startType);
        
        std::vector<std::pair<std::string, std::string>> virtualOutputsVector;
        std::vector<std::pair<std::string, std::string>> virtualInputsVector;
        std::string virtualOutputNodeName = "";
        std::string virtualInputNodeName = "";
        if(virtualNodesGranularity == "type"){
            virtualOutputNodeName = "v-out:" + endType;
            virtualInputNodeName = "v-in:" + startType;
        } else {
            virtualOutputNodeName = "v-out:" + endType + "_" + endNodeName;
            virtualInputNodeName = "v-in:" + startType + "_" + startNodeName;

        }
        
        // mapped virtualOutputs and mapped virtualInputs are the same in the sizes and logic, but have different names
        // mapped virtual outputs have the format (sourceNode, v-out:tTarget<_targetNode>)
        if(virtualNodesGranularity == "typeAndNode"){
            if(typesPairMappedVirtualOutputsVectors.contains(keyTypes)){
                if(!vectorContains(typesPairMappedVirtualOutputsVectors[keyTypes],std::make_pair(startNodeName, virtualOutputNodeName))){
                    typesPairMappedVirtualOutputsVectors[keyTypes].push_back(std::make_pair(startNodeName, virtualOutputNodeName));
                }
            } else {
                typesPairMappedVirtualOutputsVectors[keyTypes] = std::vector<std::pair<std::string, std::string>>();
                typesPairMappedVirtualOutputsVectors[keyTypes].push_back(std::make_pair(startNodeName, virtualOutputNodeName));
            }
            
            // mapped virtual input have the format (v-in:tSource<_sourceNode>, targetNode)
            if(typesPairMappedVirtualInputsVectors.contains(keyTypes)){
                if(!vectorContains(typesPairMappedVirtualInputsVectors[keyTypes],std::make_pair(virtualInputNodeName, endNodeName))){
                    typesPairMappedVirtualInputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, endNodeName));
                }
            } else {
                typesPairMappedVirtualInputsVectors[keyTypes] = std::vector<std::pair<std::string, std::string>>();
                typesPairMappedVirtualInputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, endNodeName));
            }
        } else if (virtualNodesGranularity == "type"){
            if(typesPairMappedVirtualOutputsVectors.contains(keyTypes)){
                if(!vectorContains(typesPairMappedVirtualOutputsVectors[keyTypes],std::make_pair(virtualInputNodeName, virtualOutputNodeName))){
                    typesPairMappedVirtualOutputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, virtualOutputNodeName));
                }
            } else {
                typesPairMappedVirtualOutputsVectors[keyTypes] = std::vector<std::pair<std::string, std::string>>();
                typesPairMappedVirtualOutputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, virtualOutputNodeName));
            }
            
            // mapped virtual input have the format (v-in:tSource<_sourceNode>, targetNode)
            if(typesPairMappedVirtualInputsVectors.contains(keyTypes)){
                if(!vectorContains(typesPairMappedVirtualInputsVectors[keyTypes],std::make_pair(virtualInputNodeName, virtualOutputNodeName))){
                    typesPairMappedVirtualInputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, virtualOutputNodeName));
                }
            } else {
                typesPairMappedVirtualInputsVectors[keyTypes] = std::vector<std::pair<std::string, std::string>>();
                typesPairMappedVirtualInputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, virtualOutputNodeName));
            }
        }
    }

    // populate the ranks map
    for(int rankTarget = 0; rankTarget < numProcesses; rankTarget++){
        int targetWorkload;
        if(rankTarget == numProcesses - 1){
            targetWorkload = types.size() - (numProcesses - 1)*workloadPerProcess;
        } else {
            targetWorkload = workloadPerProcess;
        }
        for(int sourceTarget = 0; sourceTarget < numProcesses; sourceTarget++){
            int sourceWorkload;
            if(sourceTarget == numProcesses - 1){
                sourceWorkload = types.size() - (numProcesses - 1)*workloadPerProcess;
            } else {
                sourceWorkload = workloadPerProcess;
            }
            for(int sourceIndexLocal = 0; sourceIndexLocal < sourceWorkload; sourceIndexLocal++){
                int sourceIndexGlobal = sourceIndexLocal + sourceTarget*workloadPerProcess;
                std::string sourceType = types[sourceIndexGlobal];
                for(int targetIndexLocal = 0; targetIndexLocal < targetWorkload; targetIndexLocal++){
                    int targetIndexGlobal = targetIndexLocal + rankTarget*workloadPerProcess;
                    std::string targetType = types[targetIndexGlobal];
                    std::pair<std::string, std::string> keyTypes = std::make_pair(sourceType, targetType);
                    std::pair<int, int> keyRanks = std::make_pair(sourceTarget, rankTarget);
                    // mapped virtual nodes have the format (v-out:tTarget<_targetNode> , v-in:tSource<_sourceNode>)
                    if(typesPairMappedVirtualOutputsVectors.contains(keyTypes)){
                        for(uint index = 0; index < typesPairMappedVirtualOutputsVectors[keyTypes].size(); index++){
                            std::pair<std::string, std::string> virtualOutput = typesPairMappedVirtualOutputsVectors[keyTypes][index];
                            std::pair<std::string, std::string> virtualInput = typesPairMappedVirtualInputsVectors[keyTypes][index];
                            std::pair<std::string, std::string> virtualNode = std::make_pair(virtualOutput.second, virtualInput.first);
                            
                            ranksPairMappedVirtualNodesVectors[keyRanks].push_back(virtualNode);
                        }
                    }
                    
                }
            }
        }
        
    }
    

    // setting propagation model in this moment since in the case of the original model, the pseudoinverse should be computed for the augmented pathway
    std::function<double(double)> propagationScalingFunction = [](double time)->double{return 1;};
    if(vm.count("propagationModel")){
        logger << "[LOG] propagation model was set to "
    << vm["propagationModel"].as<std::string>() << ".\n";
        std::string propagationModelName = vm["propagationModel"].as<std::string>();
        if(propagationModelName == "default"){
            logger << "[LOG] propagation model set to default (pseudoinverse propagation)\n";
            for(int i = 0; i < finalWorkload ;i++ ){
                typeComputations[i]->setPropagationModel(new PropagationModelOriginal(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction));
            }
            //nothing to do, default propagation scaling function is the identity
        } else if (propagationModelName == "scaled"){
            if (vm.count("propagationModelParameters")) {
                logger << "[LOG] propagation model parameters were declared to be "
            << vm["propagationModelParameters"].as<std::vector<double>>()[0] << ".\n";
                std::vector<double> propagationModelParameters = vm["propagationModelParameters"].as<std::vector<double>>();
                if(propagationModelParameters.size() == 1){
                    propagationScalingFunction = [propagationModelParameters](double time)->double{return propagationModelParameters[0];};
                    for(int i = 0; i < finalWorkload ;i++ ){
                        PropagationModel* tmpPropagationModel = new PropagationModelOriginal(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction);
                        typeComputations[i]->setPropagationModel(tmpPropagationModel);
                    }
                } else {
                    std::cerr << "[ERROR] propagation model parameters for scaled propagation must be one parameter: aborting"<<std::endl;
                    return 1;
                }
            } else {
                std::cerr << "[ERROR] propagation model parameters for scaled propagation was not set: setting to default 1 costant"<<std::endl;
                for(int i = 0; i < finalWorkload ;i++ ){
                    PropagationModel* tmpPropagationModel = new PropagationModelOriginal(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction);
                    typeComputations[i]->setPropagationModel(tmpPropagationModel);
                }
                //nothing to do, default propagation scaling function is the identity
            }
        } else if (propagationModelName == "neighbors"){
            if (vm.count("propagationModelParameters")) {
                logger << "[LOG] propagation model parameters were declared to be "
            << vm["propagationModelParameters"].as<std::vector<double>>()[0] << ".\n";
                std::vector<double> propagationModelParameters = vm["propagationModelParameters"].as<std::vector<double>>();
                if(propagationModelParameters.size() == 1){
                    propagationScalingFunction = [propagationModelParameters](double time)->double{return propagationModelParameters[0];};
                    for(int i = 0; i < finalWorkload;i++ ){
                        PropagationModel* tmpPropagationModel = new PropagationModelNeighbors(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction);
                        typeComputations[i]->setPropagationModel(tmpPropagationModel);
                    }
                } else {
                    std::cerr << "[ERROR] propagation model parameters for scaled propagation must be one parameter: aborting"<<std::endl;
                    return 1;
                }
            } else {
                std::cerr << "[ERROR] propagation model parameters for scaled propagation was not set: setting to default 1 costant"<<std::endl;
                for(int i = 0; i < finalWorkload;i++ ){
                    PropagationModel* tmpPropagationModel = new PropagationModelNeighbors(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction);
                    typeComputations[i]->setPropagationModel(tmpPropagationModel);
                }
                //nothing to do, default propagation scaling function is the identity
            }
        } else if(propagationModelName == "customScaling"){ 
            if(vm.count("propagationModelParameters")){
                logger << "[LOG] propagation model parameters were declared to be "
                << vm["propagationModelParameters"].as<std::vector<double>>()[0] << ", these parameters are not used since the propagation scaling function was set to custom.\n";
            }
            propagationScalingFunction = getPropagationScalingFunction();
            for(int i = 0; i < finalWorkload;i++ ){
                PropagationModel* tmpPropagationModel = new PropagationModelOriginal(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction);
                typeComputations[i]->setPropagationModel(tmpPropagationModel);
            }
        } else if(propagationModelName == "customScalingNeighbors"){ 
            if(vm.count("propagationModelParameters")){
                logger << "[LOG] propagation model parameters were declared to be "
                << vm["propagationModelParameters"].as<std::vector<double>>()[0] << ", these parameters are not used since the propagation scaling function was set to custom.\n";
            }
            propagationScalingFunction = getPropagationScalingFunction();
            for(int i = 0; i < finalWorkload;i++ ){
                PropagationModel* tmpPropagationModel = new PropagationModelNeighbors(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction);
                typeComputations[i]->setPropagationModel(tmpPropagationModel);
            }
            
        } else if(propagationModelName == "customPropagation"){
            if(vm.count("propagationModelParameters")){
                logger << "[LOG] propagation model parameters were declared to be "
                << vm["propagationModelParameters"].as<std::vector<double>>()[0] << ", these parameters are not used since the propagation scaling function was set to custom.\n";
            }
            propagationScalingFunction = getPropagationScalingFunction();
            for(int i = 0; i < finalWorkload;i++ ){
                PropagationModel* tmpPropagationModel = new PropagationModelCustom(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction);
                typeComputations[i]->setPropagationModel(tmpPropagationModel);
            }
        
        } else {
            std::cerr << "[ERROR] propagation model is not any of the types. propagation model scale functions available are default, scaled, neighbors and custom \n";
            return 1;
        }
    } else {
        logger << "[LOG] propagation model was not set. set to default (none)\n";
        for(int i = 0; i < finalWorkload;i++ ){
            PropagationModel* tmpPropagationModel = new PropagationModelOriginal(typeComputations[i]->getAugmentedGraph(),propagationScalingFunction);
            typeComputations[i]->setPropagationModel(tmpPropagationModel);
        }
    }

    // virtual inputs and virtual outputs buffers for the MPI communication
    // buffer for the virtual outputs from the other processes(virtual inputs), the maximum size is the power of 2 of the workload per process(since every type will send values to every other type)
    std::vector<double*> rankVirtualInputsBuffer;
    std::vector<uint> rankVirtualInputsSizes = std::vector<uint>(numProcesses,0);
    if(virtualNodesGranularity == "type"){
        for(int i = 0; i < numProcesses; i++){
            int currentWorkload;
            if(i == (numProcesses-1)){
                currentWorkload = types.size() - (i*workloadPerProcess);
            } else {
                currentWorkload = workloadPerProcess;
            }
            rankVirtualInputsSizes[i] = finalWorkload * currentWorkload;
            rankVirtualInputsBuffer.push_back(new double[finalWorkload * currentWorkload]);   // the array contains all the virtual outputs for the process types
        }
    } else if (virtualNodesGranularity == "typeAndNode"){
        // allocate the buffer for the virtual outputs for the combination of types and nodes
        for(int sourceRank = 0; sourceRank < numProcesses; sourceRank++){
            std::pair<int, int> keyRanks = std::make_pair(sourceRank, rank);
            if(ranksPairMappedVirtualNodesVectors.contains(keyRanks)){
                rankVirtualInputsSizes[sourceRank] = ranksPairMappedVirtualNodesVectors[keyRanks].size();
            } else {
                rankVirtualInputsSizes[sourceRank] = 0;
            }
            // allocate the array for the virtual inputs directed to the local rank types
            rankVirtualInputsBuffer.push_back(new double[rankVirtualInputsSizes[sourceRank]]);
        }
            
    }

    // buffer for virtual inputs to other processes(virtual outputs)
    std::vector<double*> virtualOutputs;
        std::vector<uint> rankVirtualOutputsSizes = std::vector<uint>(numProcesses,0);
        // different granularity for the virtual nodes means different ways of building the virtual outputs arrays and sizes
        // TODO generalize, difficult though since the "type" granularity has a fixed number of spots for each type, while the "typeAndNode" granularity has a variable number of spots for each type
        if(virtualNodesGranularity == "type"){ //classical way of building the virtual outputs arrays, one array for each type representing virtual nodes for each type
            for(int i = 0; i < numProcesses; i++){
                int currentWorkload;
                if(i == (numProcesses-1)){
                    currentWorkload = types.size() - (i*workloadPerProcess);
                } else {
                    currentWorkload = workloadPerProcess;
                }
                rankVirtualOutputsSizes[i] = currentWorkload * finalWorkload;
                virtualOutputs.push_back(new double[rankVirtualOutputsSizes[i]]);   // the array contains all the virtual outputs for the process types
            }
        } else if (virtualNodesGranularity == "typeAndNode"){ // finer granularity, one array for each type and node representing virtual nodes for each type and node (as a couple)
            
            //allocate the virtual outputs arrays
            for(int targetRank = 0; targetRank < numProcesses; targetRank++){
                std::pair<int, int> keyRanks = std::make_pair(rank, targetRank);
                if(ranksPairMappedVirtualNodesVectors.contains(keyRanks)){
                    rankVirtualOutputsSizes[targetRank] = ranksPairMappedVirtualNodesVectors[keyRanks].size();
                } else {
                    rankVirtualOutputsSizes[targetRank] = 0;
                }

                // allocate the array for the virtual outputs directed to the target rank types
                virtualOutputs.push_back(new double[rankVirtualOutputsSizes[targetRank]]);                        
            }
        } else {
            std::cerr << "[ERROR] virtual nodes granularity is not any of the types. virtual nodes granularity available are type and typeAndNode \n";
            return 1;
        } 

    for(int iterationInterType = 0; iterationInterType < intertypeIterations; iterationInterType++){
        for(int iterationIntraType = 0; iterationIntraType < intratypeIterations; iterationIntraType++){
            
            // computation of perturbation
            #pragma omp parallel for
            for(int i = 0; i < finalWorkload; i++){
                logger << "[LOG] computation of perturbation for iteration intertype-intratype ("+ std::to_string(iterationInterType) + "<->"+ std::to_string(iterationIntraType) + ") for type (" + types[i+startIdx]<<std::endl; 
                // TODO use stateful scaling function to consider previous times
                if (saturation) {
                    if(vm.count("saturationTerm") == 0){
                        // std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced2((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = true);
                        //std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced3((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = true, std::vector<double>(), std::vector<double>(), propagationScalingFunction);
                        std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced4((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = true);
                    } else if (vm.count("saturationTerm") >= 1) {
                        double saturationTerm = vm["saturationTerm"].as<double>();
                        std::vector<double> saturationVector = std::vector<double>(graphsNodes[invertedTypesIndexes[i]].size(),saturationTerm);
                        // std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced2((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = true, saturationVector);
                        //std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced3((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = true, saturationVector, std::vector<double>(), propagationScalingFunction); 
                        std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced4((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = true, saturationVector);
                    }
                } else{
                    // std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced2((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = false);
                    //std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced3((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = false, std::vector<double>(), std::vector<double>(), propagationScalingFunction);
                    std::vector<double> outputValues = typeComputations[i]->computeAugmentedPerturbationEnhanced4((iterationInterType*intratypeIterations + iterationIntraType)*timestep, saturation = false);
                }
            }



            //save output values
            for(int i = 0; i < finalWorkload; i++){
                std::vector<std::string> nodeNames = typeComputations[i]->getAugmentedGraph()->getNodeNames();
                //TODO change how to save files to get more information about intratype and intertype iterations
                //logger << "saving output values for iteration intertype-intratype ("+ std::to_string(iterationInterType) + "<->"+ std::to_string(iterationIntraType) + ") for type (" + types[i+startIdx] << ") in process " << rank <<std::endl;
                saveNodeValuesWithTime(outputFoldername, iterationInterType*intratypeIterations, iterationIntraType, types[i+startIdx], typeComputations[i]->getOutputAugmented(), nodeNames,ensembleGeneNames, nodesDescriptionFilename, timestep);
            }

            //update input
            for(int i = 0; i < finalWorkload; i++){
                //If conservation of the initial values is required, the input is first updated with the initial norm value
                if (conservateInitialNorm) {
                    int index = indexMapGraphTypesToValuesTypes[i+startIdx];
                    std::vector<double> inputInitial = inputInitials[index];
                    double initialNorm = vectorNorm(inputInitial);
                    double outputNorm = vectorNorm(typeComputations[i]->getOutputAugmented());
                    double normRatio = initialNorm/outputNorm;
                    std::vector<double> newInput = vectorScalarMultiplication(typeComputations[i]->getOutputAugmented(),normRatio);
                    logger << "[LOG] update input with conservation of the initial perturbation for iteration intertype-intratype ("+ std::to_string(iterationInterType) + "<->"+ std::to_string(iterationIntraType) + ") for type (" + types[i+startIdx]<<std::endl;
                    typeComputations[i]->updateInput(newInput,true);
                } else {
                    logger << "[LOG] update input for iteration intertype-intratype ("+ std::to_string(iterationInterType) + "<->"+ std::to_string(iterationIntraType) + ") for type (" + types[i+startIdx]<<std::endl;
                    typeComputations[i]->updateInput(std::vector<double>(),true);
                }
                
            }
        }

        // send virtual outputs to the other processes, the vector contains the virtual outputs for each type as an array
        // for every type, send the virtual outputs to the other processes, all in the same array (this array will be decomposed on the target)
        // build the array


        if(virtualNodesGranularity == "type"){ //classical way of building the virtual outputs arrays, one array for each type representing virtual nodes for each type
            // fill the arrays
            for(int i = 0; i < SizeToInt(types.size()); i++){
                int targetRank = typeToRank[types[i]];
                int targetWorkload;
                if(targetRank == (numProcesses-1)){
                    targetWorkload = types.size() - (targetRank*workloadPerProcess);
                } else {
                    targetWorkload = workloadPerProcess;
                }
                    int targetPosition = i - targetRank * workloadPerProcess;
                    for(int j = 0; j < finalWorkload; j++ ){
                        int virtualOutputPosition = targetPosition + j * targetWorkload;
                        
                        virtualOutputs[targetRank][virtualOutputPosition] = typeComputations[j]->getVirtualOutputForType(types[i]);
                    }
                
            }
        } else if (virtualNodesGranularity == "typeAndNode"){ // finer granularity, one array for each type and node representing virtual nodes for each type and node (as a couple)    
            // fill the arrays
            for (int targetRank = 0; targetRank < numProcesses; targetRank++){
                std::pair<int, int> keyRanks = std::make_pair(rank, targetRank);
                if(ranksPairMappedVirtualNodesVectors.contains(keyRanks)){
                    for(uint i = 0; i < ranksPairMappedVirtualNodesVectors[keyRanks].size(); i++){
                        std::pair<std::string, std::string> virtualNode = ranksPairMappedVirtualNodesVectors[keyRanks][i];
                        std::vector<std::string> virtualInputNodeSplit = splitVirtualNodeStringIntoVector(virtualNode.second);
                        if(virtualInputNodeSplit.size()==3){
                            std::string sourceType = virtualInputNodeSplit[1];
                            std::string sourceNode = virtualInputNodeSplit[2];
                            int sourceTypePosition = -1;
                            for(int j = 0; j < finalWorkload; j++){
                                if(types[j+startIdx] == sourceType){
                                    sourceTypePosition = j;
                                    break;
                                }
                            }
                            if(sourceTypePosition == -1){
                                std::cerr << "[ERROR] source type not found in the types vector: aborting" << std::endl;
                                return 1;
                            }
                            double virtualOutputValue = typeComputations[sourceTypePosition]->getOutputNodeValue(virtualNode.first);
                            virtualOutputs[targetRank][i] = virtualOutputValue;
                        } else {
                            std::cerr << "[ERROR] virtual node string is not in the correct format: aborting" << std::endl;
                            return 1;
                        }
                    }
                }
                        
            }
        } else {
            std::cerr << "[ERROR] virtual nodes granularity is not any of the types. virtual nodes granularity available are type and typeAndNode \n";
            return 1;
        }


        


        // reset virtual outputs if specified, should work since virtual outputs are assigned before
        if(resetVirtualOutputs){
            for(int i = 0; i < finalWorkload; i++){
                typeComputations[i]->resetVirtualOutputs();
            }
        }



        // preliminary asynchronous receive
        
        // receive the virtual outputs from the other processes
        MPI_Request request[numProcesses];
        for(int i = 0; i < numProcesses; i++){
            int sourceRank = i;
            // receive only the virtual outputs for the types granularity (v-out for each type), or the full vectors for each pair of source type and target type
            if(virtualNodesGranularity == "typeAndNode" || virtualNodesGranularity == "type"){
                std::pair<int, int> keyRanks = std::make_pair(sourceRank, rank);
                if(ranksPairMappedVirtualNodesVectors.contains(keyRanks)){
                    MPI_Irecv(rankVirtualInputsBuffer[sourceRank], rankVirtualInputsSizes[sourceRank], MPI_DOUBLE, sourceRank, 0, MPI_COMM_WORLD, &request[sourceRank]);
                }
            } else {
                // OTHER CASES NOT IMPLEMENTED YET
            }
        }

    

        // send the virtual outputs to the other processes
        for(int targetRank = 0; targetRank < numProcesses; targetRank++){
            //sending virtual outputs to target cell
            // int targetStartIdx = j * workloadPerProcess;
            // int targetEndIdx = (j == numProcesses - 1) ? types.size() : (j + 1) * workloadPerProcess;
            // logger << "[LOG] sending virtual output from type " << types[startIdx] << " to type " << types[endIdx-1] << " from process " << rank << " to process " << j << " from type " << types[targetStartIdx] << " to type " << types[targetEndIdx-1] << std::endl;
            
            //synchronized communication will lead to deadlocks with this type of implementation
            if(virtualNodesGranularity == "typeAndNode" || virtualNodesGranularity == "type"){
                //send the subvectors of the virtual outputs for the combination of types and nodes
                std::pair<int, int> keyRanks = std::make_pair(rank, targetRank);
                if(ranksPairMappedVirtualNodesVectors.contains(keyRanks)){
                    try
                    {
                        MPI_Send(virtualOutputs.at(targetRank), rankVirtualOutputsSizes[targetRank], MPI_DOUBLE, targetRank, 0, MPI_COMM_WORLD);
                    }
                    catch(const std::exception& e)
                    {
                        std::cerr << e.what() << std::endl;
                        std::cerr << "[ERROR] error in sending virtual outputs from process " << rank << " to process " << targetRank << std::endl;
                        return 1;
                    }
                    logger << "[LOG] sent virtual outputs from process " << rank << " to process " << targetRank  << std::endl;
                }
            } else {
                // send only the virtual outputs for the types granularity (v-out for each type
            }
        }

        // receive outputs from the other processes and update the input
        for(int sourceRank = 0; sourceRank < numProcesses; sourceRank++){
            std::pair<int, int> ranksPair = std::make_pair(sourceRank, rank);
            if(ranksPairMappedVirtualNodesVectors.contains(ranksPair)){
                logger << "[LOG] receiving virtuals outputs from process " << sourceRank << " to process " << rank << std::endl;
                try{
                    MPI_Wait(&request[sourceRank], MPI_STATUS_IGNORE);
                } catch(const std::exception& e){
                    std::cerr << e.what() << std::endl;
                    std::cerr << "[ERROR] error in waiting for virtual outputs from process " << sourceRank << " to process " << rank << std::endl;
                    return 1;
                }
                logger << "[LOG] received virtual outputs from process " << sourceRank << " to process " << rank << std::endl;
            }
            // source workload and virtual outputs decomposition on the target

            int sourceWorkload;
            if(sourceRank == (numProcesses-1)){
                sourceWorkload = types.size() - (sourceRank*workloadPerProcess);
            } else {
                sourceWorkload = workloadPerProcess;
            }
            if(virtualNodesGranularity == "type"){
                for(int isource = 0; isource < sourceWorkload; isource++){
                    for(int ilocal = 0; ilocal < finalWorkload; ilocal++){
                        int virtualInputPosition = ilocal + isource * finalWorkload;
                        int localTypePosition = ilocal + startIdx;
                        int sourceTypePosition = isource + sourceRank*workloadPerProcess;
                        std::pair<std::string,std::string> keyTypes = std::make_pair(types[localTypePosition], types[sourceTypePosition]);
                        if(interactionBetweenTypesMap[keyTypes].contains(iterationInterType)){
                            // logger << "[TEST] contact times for types " << types[localTypePosition] << " and " << types[sourceTypePosition] << " are ";
                            // for(auto time: interactionBetweenTypesMap[keyTypes]){
                            //     logger << time << ", ";
                            // } 
                            // logger << std::endl;
                            if(localTypePosition==sourceTypePosition){
                                if(sameTypeCommunication) typeComputations[ilocal]->setInputVinForType(types[sourceTypePosition], rankVirtualInputsBuffer[sourceRank][virtualInputPosition]);
                            } else {
                                typeComputations[ilocal]->setInputVinForType(types[sourceTypePosition], rankVirtualInputsBuffer[sourceRank][virtualInputPosition]);
                            }
                        }
                    }
                }
            } else if (virtualNodesGranularity == "node"){
                logger << "[ERROR] virtual nodes granularity is not supported yet: aborting"<<std::endl;
                return 1;
            } else if (virtualNodesGranularity == "typeAndNode"){
                // logic of reading the subvectors of the virtual inputs

                // use ranksPairMappedVirtualNodesVectors to get all the information needed to update the input
                int targetRank = rank;
                std::pair<int,int> ranksPair = std::make_pair(sourceRank,targetRank);

                if(ranksPairMappedVirtualNodesVectors.contains(ranksPair)){

                    for(uint i = 0; i < ranksPairMappedVirtualNodesVectors[ranksPair].size(); i++){
                        
                        std::pair<std::string, std::string> virtualNodes = ranksPairMappedVirtualNodesVectors[ranksPair][i];
                        std::string virtualOutputNodeName = virtualNodes.first;
                        std::string virtualInputNodeName = virtualNodes.second;
                        std::vector<std::string> virtualOutputNodeNameSplitted = splitVirtualNodeStringIntoVector(virtualOutputNodeName);
                        if(virtualOutputNodeNameSplitted.size()!=3) throw std::runtime_error("main:: virtual output node name is not in the correct format: " + virtualOutputNodeName);
                        std::string targetType = virtualOutputNodeNameSplitted[1];
                        std::string targetNodeName = virtualOutputNodeNameSplitted[2];

                        std::vector<std::string> virtualInputNodeNameSplitted = splitVirtualNodeStringIntoVector(virtualInputNodeName);
                        if(virtualInputNodeNameSplitted.size()!=3) throw std::runtime_error("main:: virtual input node name is not in the correct format: " + virtualInputNodeName);
                        std::string sourceType = virtualInputNodeNameSplitted[1];
                        std::string sourceNodeName = virtualInputNodeNameSplitted[2];

                        // get local target index for typeComputation
                        int targetTypeIndex = -1;
                        for(int i = 0; i < finalWorkload; i++){
                            if(types[i+startIdx] == targetType){
                                targetTypeIndex = i;
                                break;
                            }
                        }

                        if(targetTypeIndex == -1) throw std::runtime_error("main:: target type index not found for type: " + targetType);
                        std::tuple<std::string, std::string, std::string, std::string> interactionKey = std::make_tuple(sourceNodeName, targetNodeName, sourceType, targetType);
                        if(interactionBetweenTypesFinerMap.contains(interactionKey)){
                            if(interactionBetweenTypesFinerMap[interactionKey].contains(iterationInterType)){
                                double newValue = rankVirtualInputsBuffer[sourceRank][i];
                                try{
                                    if(sourceType == targetType){
                                        if(sameTypeCommunication) typeComputations[targetTypeIndex ]->setInputNodeValue(virtualInputNodeName, newValue);
                                    } else {
                                        typeComputations[targetTypeIndex]->setInputNodeValue(virtualInputNodeName, newValue);
                                    }
                                } catch(const std::exception& e){
                                    std::cerr << e.what() << std::endl;
                                    std::cerr << "[ERROR] error in setting input for virtual nodes from process "<< sourceRank<< " to process "<< targetRank << " for virtual nodes: " << virtualInputNodeName << " and " << virtualOutputNodeName << std::endl;
                                    return 1;
                                }
                            } else {
                                // // TESTING
                                // std::cout << "[DEBUG] rank: " << rank << " interaction between nodes " << sourceNodeName << " and " << targetNodeName << " for types " << sourceType << " and " << targetType << " have no contact time for inter iteration "<< iterationInterType << std::endl;
                                // // TESTING
                            }
                        } else {
                            std::cerr << "[ERROR] interaction between nodes " << sourceNodeName << " and " << targetNodeName << " for types " << sourceType << " and " << targetType << " is not present in the interactionBetweenTypesFinerMap" << std::endl;
                            std::cerr << "[ERROR] aborting" << std::endl;
                            return 1;
                        }

                        
                    }
                }
            }
        }

        // // TESTING
        // // printing nodes values(inputAugmented) for augmented graph 1
        // for(int i = 0; i< finalWorkload; i++){
        //     if(types[startIdx + i] == "1"){
        //         std::vector<std::string> nodeNames = typeComputations[i]->getAugmentedGraph()->getNodeNames();
        //         std::cout << "printing values for type 1 after receiving them via MPI"<<std::endl;
        //         for(uint nodeIndex = 0; nodeIndex < nodeNames.size(); nodeIndex++){
        //             std::cout << "(" << nodeNames[nodeIndex] << " = " << typeComputations[i]->getInputAugmented()[nodeIndex] << "), ";
        //         }
        //     }
        // } 
        // std::cout << std::endl;
        // // TESTING
    }
    // delete the virtual outputs vector of arrays
    for(uint i = 0; i < virtualOutputs.size(); i++){
        if (virtualOutputs.at(i) != nullptr) delete[] virtualOutputs.at(i);
    }
    // delete the virtual inputs vector of arrays
    for(uint i = 0; i < rankVirtualInputsBuffer.size(); i++){
        if (rankVirtualInputsBuffer.at(i) != nullptr) delete[] rankVirtualInputsBuffer.at(i);
    }

    MPI_Finalize();
    // take ending time after the computation
    auto end = std::chrono::steady_clock::now();
    if(rank == 0){
        if(vm.count("savePerformance")){
            std::ofstream performanceFile;
            int numberProcesses = numProcesses;
            int numberTypes = types.size();
            int numberIterations = intratypeIterations * intertypeIterations;
            performanceFile.open (performanceFilename, std::ios::out | std::ios::app);
            if (performanceFile.fail())
                throw std::ios_base::failure(std::strerror(errno));
            //performanceFile << "inputFolderGraphs\t" << "numberProcesses" << "\t" << "numberTypes" << "\t" << "numberIterations" << "\t" << "time" << std::endl;
            performanceFile << graphsFilesFolder << "\t" << numberProcesses << "\t" << numberTypes << "\t" << numberIterations << "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
            performanceFile.close();
        }
    }
    return 0;
}
