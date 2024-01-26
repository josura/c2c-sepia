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
    // starting data for each process
    
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
    // // TESTING
    // std::cout << "rank: " << rank << " startIdx: " << startIdx << " endIdx: " << endIdx << " finalWorkload: " << finalWorkload << std::endl;
    // // END TESTING


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

    // // TESTING
    // std::cout << "rank: " << rank << " typeToRank: " << std::endl;
    // for (auto const& x : typeToRank)
    // {
    //     std::cout << x.first  // string (key)
    //               << ':'
    //               << x.second // string's value 
    //               << std::endl ;
    // }
    // // END TESTING


    //use the number of types for workload to allocate an array of pointers to contain the graph for each type
    WeightedEdgeGraph **graphs = new WeightedEdgeGraph*[finalWorkload];
    std::vector<std::vector<std::string>> graphsNodes;
    std::vector<std::vector<std::string>> graphsNodesAll; // used only initially to read the values, contains all types
    std::unordered_map<std::string, std::vector<std::string>> typeToNodeNamesMap; // map from all types to the node names, not only the ones in the workload
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
            typeComputations[indexComputation]->augmentGraphNoComputeInverse(types,std::vector<std::pair<std::string,std::string>>(),std::vector<double>(), true);
        }
        typesIndexes[i] = indexComputation;
        invertedTypesIndexes[indexComputation] = i;
        indexComputation++;

    }

    // read the type interactions
    std::vector<std::vector<std::string>> typeToNodeNames = std::vector<std::vector<std::string>>(finalWorkload,std::vector<std::string>());
    
    for(int i = 0; i < finalWorkload;i++ ){
        typeToNodeNames[i] = typeComputations[i]->getAugmentedGraph()->getNodeNames();    
    }
    auto allFilesInteraction = get_all(typesInteractionFoldername,".tsv");
    // define the map for the type interactions, an hash function should be defined for the pair of strings used as the identifier of the interaction
    // TODO substitute with another class that represents granularity and returns the interactions between the types or the pairs of types+node, along the lists of contact times and virtual nodes (just a superclass that is extended by the two classes for different granularity)
    std::unordered_map<std::pair<std::string, std::string>, std::unordered_set<int>, hash_pair_strings> interactionBetweenTypesMap;
    std::unordered_map<std::tuple<std::string, std::string, std::string, std::string>, std::unordered_set<int>, hash_quadruple_strings> interactionBetweenTypesFinerMap;
    for(auto typeInteractionFilename = allFilesInteraction.cbegin() ; typeInteractionFilename != allFilesInteraction.cend() ; typeInteractionFilename++){
        std::pair<std::map<std::string,std::vector<std::tuple<std::string,std::string,double>>>,std::vector<std::tuple<std::string, std::string, std::string, std::string, std::unordered_set<int>,double>>> typeInteractionsEdges;
        if (subtypes.size() == 0) {
            // TODO add different contact times inside the network (quite difficult since the structure of the graphs is static)
            // SOLUTION: granularity
            typeInteractionsEdges  = interactionContactsFileToEdgesListAndNodesByName(*typeInteractionFilename, types, intertypeIterations, ensembleGeneNames, virtualNodesGranularity, typeToNodeNamesMap);
        } else {
            typeInteractionsEdges = interactionContactsFileToEdgesListAndNodesByName(*typeInteractionFilename, subtypes, intertypeIterations, ensembleGeneNames, virtualNodesGranularity, typeToNodeNamesMap);
        }
        #pragma omp parallel for
        for (int i = 0; i < finalWorkload;i++) {
            if(typeInteractionsEdges.first.contains(types[i+startIdx]) && typesIndexes[i] != -1){
                // granularity is already considered in the function that reads from the file previously called
                typeComputations[typesIndexes[i]]->addEdges(typeInteractionsEdges.first[types[i+startIdx]], undirectedTypeEdges, false); // no inverse computation since it is done in the propagation model
            }
        }
        for(auto edge = typeInteractionsEdges.second.cbegin() ; edge != typeInteractionsEdges.second.cend(); edge++ ){
            // first two types are the nodes in the two networks/types ,types are the third and fourth element of the tuple, while the fifth is the set of contact times
            // startNodeName, endNodeName, startType, endType, contactTimes, weight
            std::pair<std::string,std::string> keyTypes = std::make_pair(std::get<2> (*edge), std::get<3> (*edge));
            std::tuple<std::string,std::string,std::string,std::string> keyTypesFiner = std::make_tuple(std::get<0> (*edge), std::get<1> (*edge), std::get<2> (*edge), std::get<3> (*edge));
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

    // create a map that maps couples of strings (source type and target type) to a vector of pairs of strings, representing how the virtual outputs are mapped in the array passed to MPI send 
    std::unordered_map<std::pair<std::string, std::string>, std::vector<std::pair<std::string, std::string>>,hash_pair_strings> mappedVirtualOutputsVectors;
    // create a map that maps couples of strings (source type and target type) to a vector of pairs of strings, representing how the virtual inputs are mapped in the array passed to MPI send
    std::unordered_map<std::pair<std::string, std::string>, std::vector<std::pair<std::string, std::string>>,hash_pair_strings> mappedVirtualInputsVectors;

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
            if(mappedVirtualOutputsVectors.contains(keyTypes)){
                if(!vectorContains(mappedVirtualOutputsVectors[keyTypes],std::make_pair(startNodeName, virtualOutputNodeName))){
                    mappedVirtualOutputsVectors[keyTypes].push_back(std::make_pair(startNodeName, virtualOutputNodeName));
                }
            } else {
                mappedVirtualOutputsVectors[keyTypes] = std::vector<std::pair<std::string, std::string>>();
                mappedVirtualOutputsVectors[keyTypes].push_back(std::make_pair(startNodeName, virtualOutputNodeName));
            }
            
            // mapped virtual input have the format (v-in:tSource<_sourceNode>, targetNode)
            if(mappedVirtualInputsVectors.contains(keyTypes)){
                if(!vectorContains(mappedVirtualInputsVectors[keyTypes],std::make_pair(virtualInputNodeName, endNodeName))){
                    mappedVirtualInputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, endNodeName));
                }
            } else {
                mappedVirtualInputsVectors[keyTypes] = std::vector<std::pair<std::string, std::string>>();
                mappedVirtualInputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, endNodeName));
            }
        } else if (virtualNodesGranularity == "type"){
            if(mappedVirtualOutputsVectors.contains(keyTypes)){
                if(!vectorContains(mappedVirtualOutputsVectors[keyTypes],std::make_pair(virtualInputNodeName, virtualOutputNodeName))){
                    mappedVirtualOutputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, virtualOutputNodeName));
                }
            } else {
                mappedVirtualOutputsVectors[keyTypes] = std::vector<std::pair<std::string, std::string>>();
                mappedVirtualOutputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, virtualOutputNodeName));
            }
            
            // mapped virtual input have the format (v-in:tSource<_sourceNode>, targetNode)
            if(mappedVirtualInputsVectors.contains(keyTypes)){
                if(!vectorContains(mappedVirtualInputsVectors[keyTypes],std::make_pair(virtualInputNodeName, virtualOutputNodeName))){
                    mappedVirtualInputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, virtualOutputNodeName));
                }
            } else {
                mappedVirtualInputsVectors[keyTypes] = std::vector<std::pair<std::string, std::string>>();
                mappedVirtualInputsVectors[keyTypes].push_back(std::make_pair(virtualInputNodeName, virtualOutputNodeName));
            }
        }
    }
    
    // setting propagation model in this moment since in the case of the original model, the pseudoinverse should be computed for the augmented pathway
    // TESTING
    if(rank == 0){
        std::cout << "[TESTING] rank: " << rank << " mapped virtual inputs and outputs" << std::endl;
        std::cout << "[TESTING] mapped virtual inputs size: " << mappedVirtualInputsVectors.size() << std::endl;
        for(auto interaction = mappedVirtualOutputsVectors.cbegin() ; interaction != mappedVirtualOutputsVectors.cend(); interaction++ ){
            std::cout << interaction->first.first << " " << interaction->first.second << " ";
            for(auto virtualOutput = interaction->second.cbegin() ; virtualOutput != interaction->second.cend(); virtualOutput++ ){
                std::cout << virtualOutput->first << " " << virtualOutput->second << " ";
            }
            std::cout << std::endl;
        }
    }
    // TESTING

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

    for(int iterationInterType = 0; iterationInterType < intertypeIterations; iterationInterType++){
        for(int iterationIntraType = 0; iterationIntraType < intratypeIterations; iterationIntraType++){
            
            
            // TESTING
            // print all the node values for t0
            for (int i = 0; i < finalWorkload; i++){
                if(types[i+startIdx] == "t0"){
                    logger << "[DEBUG] node values for type t0 in interIteration "<< iterationInterType << " and intra type iteration "<< iterationIntraType <<" before: " << std::endl;
                    std::vector<std::string> nodeNames = typeComputations[i]->getAugmentedGraph()->getNodeNames();
                    for (uint j = 0; j < nodeNames.size(); j++){
                        logger << nodeNames[j] << " = " << typeComputations[i]->getInputNodeValue(nodeNames[j]) << ", ";
                    }
                    logger << std::endl;
                    // logger << "[DEBUG] virtual outputs for type t0 in interIteration "<< iterationInterType <<" before: " << std::endl;
                    // for(int j = 0; j < numProcesses; j++){
                    //     logger << std::endl << "[DEBUG] to process "<< j << " : " << std::endl;
                    //     int targetWorkload;
                    //     if(j == (numProcesses-1)){
                    //         targetWorkload = types.size() - (j*workloadPerProcess);
                    //     } else {
                    //         targetWorkload = workloadPerProcess;
                    //     }
                    //     for(int k = 0; k < targetWorkload; k++){
                    //         int localTypePosition = i + startIdx;
                    //         int targetTypePosition = k + j*workloadPerProcess;
                    //         logger << "v(" << types[localTypePosition]<< "->" << types[targetTypePosition] << ")= v-out:"<<types[targetTypePosition]<< "= " <<typeComputations[i]->getOutputNodeValue("v-out:" + types[targetTypePosition]) << " = "  << typeComputations[i]->getVirtualOutputForType(types[targetTypePosition]) << ", ";
                    //     }
                    // }
                    // logger << std::endl;
                }
            }
            // TESTING
            
            // computation of perturbation
            #pragma omp parallel for
            for(int i = 0; i < finalWorkload; i++){
                std::vector<std::string> nodeNames = typeToNodeNames[i];
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


            // TESTING
            // print all the node values for t0
            for (int i = 0; i < finalWorkload; i++){
                if(types[i+startIdx] == "t0"){
                    logger << "[DEGUB] node values for type t0 in interIteration "<< iterationInterType << " and intra type iteration "<< iterationIntraType <<" after: " << std::endl;
                    std::vector<std::string> nodeNames = typeComputations[i]->getAugmentedGraph()->getNodeNames();
                    for (uint j = 0; j < nodeNames.size(); j++){
                        logger << nodeNames[j] << " = " << typeComputations[i]->getOutputNodeValue(nodeNames[j]) << ", ";
                    }
                    logger << std::endl;
                    // logger << "[DEBUG] virtual outputs for type t0 in interIteration "<< iterationInterType <<" before: " << std::endl;
                    // for(int j = 0; j < numProcesses; j++){
                    //     logger << std::endl << "[DEBUG] to process "<< j << " : " << std::endl;
                    //     int targetWorkload;
                    //     if(j == (numProcesses-1)){
                    //         targetWorkload = types.size() - (j*workloadPerProcess);
                    //     } else {
                    //         targetWorkload = workloadPerProcess;
                    //     }
                    //     for(int k = 0; k < targetWorkload; k++){
                    //         int localTypePosition = i + startIdx;
                    //         int targetTypePosition = k + j*workloadPerProcess;
                    //         logger << "v(" << types[localTypePosition]<< "->" << types[targetTypePosition] << ")= v-out:"<<types[targetTypePosition]<< "= " <<typeComputations[i]->getOutputNodeValue("v-out:" + types[targetTypePosition]) << " = "  << typeComputations[i]->getVirtualOutputForType(types[targetTypePosition]) << ", ";
                    //     }
                    // }
                    // logger << std::endl;
                }
            }
            // TESTING



            //save output values
            for(int i = 0; i < finalWorkload; i++){
                std::vector<std::string> nodeNames = typeToNodeNames[i];
                //TODO change how to save files to get more information about intratype and intertype iterations
                //logger << "saving output values for iteration intertype-intratype ("+ std::to_string(iterationInterType) + "<->"+ std::to_string(iterationIntraType) + ") for type (" + types[i+startIdx] << ") in process " << rank <<std::endl;
                saveNodeValues(outputFoldername, iterationInterType*intratypeIterations + iterationIntraType, types[i+startIdx], typeComputations[i]->getOutputAugmented(), nodeNames,ensembleGeneNames, nodesDescriptionFilename);
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

        // //TESTING
        // //printing type to rank map
        // logger << "[LOG] type to rank map: " << std::endl;
        // for (auto const& x : typeToRank)
        // {
        //     logger << x.first  // string (key)
        //     << ':'
        //     << x.second // string's value
        //     << ", " ;
        // }
        // logger << std::endl;
        // //printing start idx and end idx with rank
        // logger << "[LOG] start idx: " << startIdx << " end idx: " << endIdx << " rank: " << rank << std::endl;
        // // printing types
        // logger << "[LOG] types for rank "<< rank <<": " << std::endl;
        // for (auto const& x : types)
        // {
        //     logger << x << ", " ;
        // }
        // //TESTING

        // send virtual outputs to the other processes, the vector contains the virtual outputs for each type as an array
        // for every type, send the virtual outputs to the other processes, all in the same array (this array will be decomposed on the target)
        // build the array
        
        
        
        
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
                        // // TESTING
                        // logger << "[LOG] virtual output position: " << virtualOutputPosition << " for v(" << types[j + startIdx]<< "->" << types[i] << ")" << std::endl;
                        // // TESTING
                        virtualOutputs[targetRank][virtualOutputPosition] = typeComputations[j]->getVirtualOutputForType(types[i]);
                    }
                
            }
        } else if (virtualNodesGranularity == "typeAndNode"){ // finer granularity, one array for each type and node representing virtual nodes for each type and node (as a couple)
            
            //allocate the virtual outputs arrays
            for(int targetRank = 0; targetRank < numProcesses; targetRank++){
                for(int sourceIndexLocal = 0; sourceIndexLocal < finalWorkload; sourceIndexLocal++){
                    std::string sourceType = types[sourceIndexLocal + startIdx];
                    int targetWorkload;
                    int targetStartIdx = targetRank * workloadPerProcess;
                    if(targetRank == (numProcesses-1)){
                        targetWorkload = types.size() - (targetRank*workloadPerProcess);
                    } else {
                        targetWorkload = workloadPerProcess;
                    }
                    for(int targetIndexLocal = 0; targetIndexLocal < targetWorkload; targetIndexLocal++){
                        std::string targetType = types[targetIndexLocal + targetStartIdx];
                        // if there is at least an interaction between the two types, the size is increased
                        if(mappedVirtualOutputsVectors.contains(std::make_pair(sourceType,targetType))){
                            rankVirtualOutputsSizes[targetRank] += mappedVirtualOutputsVectors[std::make_pair(sourceType,targetType)].size();
                        }
                    }
                }

                // allocate the array for the virtual outputs directed to the target rank types
                virtualOutputs.push_back(new double[rankVirtualOutputsSizes[targetRank]]);

            }

            // fill the arrays
            for (int targetRank = 0; targetRank < numProcesses; targetRank++){
                int targetStartIdx = targetRank * workloadPerProcess;
                int targetWorkload;
                if(targetRank == (numProcesses-1)){
                    targetWorkload = types.size() - (targetRank*workloadPerProcess);
                } else {
                    targetWorkload = workloadPerProcess;
                }
                // code is generalized for both cases of type granularity and typeAndNode granularity
                int virtualOutputPosition = 0;
                for (int sourceIndexLocal = 0; sourceIndexLocal < finalWorkload; sourceIndexLocal++){
                    std::string sourceType = types[sourceIndexLocal + startIdx];
                    for (int targetIndexLocal = 0; targetIndexLocal < targetWorkload; targetIndexLocal++){
                        std::string targetType = types[targetIndexLocal + targetStartIdx];
                        // if there is at least an interaction between the two types, the size is increased
                        if(mappedVirtualOutputsVectors.contains(std::make_pair(sourceType,targetType))){
                            for(auto virtualOutputPair : mappedVirtualOutputsVectors[std::make_pair(sourceType,targetType)]){
                                // the virtual output is ordered by the sequence (t_i,t_j)(t_i,t_j+1)...(t_i,t_n)(t_i+1,t_j)(t_i+1,t_j+1)...(t_i+1,t_n)...(t_n,t_n)
                                // that is the sequence of the virtual outputs for the target type, for each source type
                                double virtualOutputValue = typeComputations[sourceIndexLocal]->getOutputNodeValue(virtualOutputPair.second);
                                virtualOutputs[targetRank][virtualOutputPosition] = virtualOutputValue;
                                virtualOutputPosition++;
                            }
                        }
                    }
                }
                        
            }
        } else {
            std::cerr << "[ERROR] virtual nodes granularity is not any of the types. virtual nodes granularity available are type and typeAndNode \n";
            return 1;
        }

        


        // reset virtual outputs if specified
        if(resetVirtualOutputs){
            for(int i = 0; i < finalWorkload; i++){
                typeComputations[i]->resetVirtualOutputs();
            }
        }


        

        // // TESTING
        // // print the virtual outputs arrays
        // logger << "[LOG] virtual outputs for process "<< rank<< "  after: " << std::endl;
        // for(int i = 0; i < numProcesses; i++){
        //     logger << std::endl << "[LOG] to process "<< i << " : " << std::endl;
        //     int targetWorkload;
        //     if(i == (numProcesses-1)){
        //         targetWorkload = types.size() - (i*workloadPerProcess);
        //     } else {
        //         targetWorkload = workloadPerProcess;
        //     }
        //     for(int j = 0; j < finalWorkload; j++){
        //         for(int k = 0; k < targetWorkload; k++){
        //             int virtualOutputPosition = j + k * finalWorkload;
        //             int localTypePosition = j + startIdx;
        //             int targetTypePosition = k + i*workloadPerProcess;
        //             logger << "v(" << types[localTypePosition]<< "->" << types[targetTypePosition] << ")=" << virtualOutputs.at(i)[virtualOutputPosition] << ", ";
        //         }
        //     }
        // }
        // logger << std::endl;
        // // TESTING

        


        // preliminary asynchronous receive
        // buffer for the virtual outputs from the other processes, the maximum size is the power of 2 of the workload per process(since every type will send values to every other type)
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
            // TODO allocate the buffer for the virtual outputs for the combination of types and nodes
            for(int sourceRank = 0; sourceRank < numProcesses; sourceRank++){
                // compute the virtual inputs sizes for the source rank types to the local rank types
                int sourceStartIdx = sourceRank * workloadPerProcess;
                int sourceWorkload;
                if(sourceRank == (numProcesses-1)){
                    sourceWorkload = types.size() - (sourceRank*workloadPerProcess);
                } else {
                    sourceWorkload = workloadPerProcess;
                }
                for(int targetIndexLocal = 0; targetIndexLocal < finalWorkload; targetIndexLocal++){
                    std::string targetType = types[targetIndexLocal + startIdx];
                    for(int sourceIndexLocal = 0; sourceIndexLocal < sourceWorkload; sourceIndexLocal++){
                        std::string sourceType = types[sourceIndexLocal + sourceStartIdx];
                        // if there is at least an interaction between the two types, the size is increased
                        if(mappedVirtualInputsVectors.contains(std::make_pair(sourceType,targetType))){
                            rankVirtualInputsSizes[sourceRank] += mappedVirtualInputsVectors[std::make_pair(sourceType,targetType)].size();
                        }
                    }
                }
                rankVirtualInputsBuffer.push_back(new double[rankVirtualInputsSizes[sourceRank]]);
            }
                
        }
        // receive the virtual outputs from the other processes
        MPI_Request request[numProcesses];
        for(int i = 0; i < numProcesses; i++){
            int sourceRank = i;
            // receive only the virtual outputs for the types granularity (v-out for each type), or the full vectors for each pair of source type and target type
            if(virtualNodesGranularity == "typeAndNode" || virtualNodesGranularity == "type"){
                MPI_Irecv(rankVirtualInputsBuffer[sourceRank], rankVirtualInputsSizes[sourceRank], MPI_DOUBLE, sourceRank, 0, MPI_COMM_WORLD, &request[sourceRank]);
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
                MPI_Send(virtualOutputs.at(targetRank), rankVirtualOutputsSizes[targetRank], MPI_DOUBLE, targetRank, 0, MPI_COMM_WORLD);
                logger << "[LOG] sent virtual outputs from process " << rank << " to process " << targetRank  << std::endl;
            } else {
                // send only the virtual outputs for the types granularity (v-out for each type
            }
        }

        // delete the virtual outputs vector of arrays
        for(int i = 0; i < numProcesses; i++){
            delete[] virtualOutputs.at(i);
        }

        // receive outputs from the other processes and update the input
        for(int sourceRank = 0; sourceRank < numProcesses; sourceRank++){
            //logger << "[LOG] receiving virtuals outputs from process " << sourceRank << " to process " << rank << std::endl;
            MPI_Wait(&request[sourceRank], MPI_STATUS_IGNORE);
            logger << "[LOG] received virtual outputs from process " << sourceRank << " to process " << rank << std::endl;
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
                        // TODO take into account granularity of virtual nodes in the future
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
                // TODO implement the logic of reading the subvectors of the virtual outputs   
                // int sourceStartIdx = sourceRank * workloadPerProcess;
                // int virtualInputPosition = 0;
                // for(int targetIndexLocal = 0; targetIndexLocal < finalWorkload; targetIndexLocal++){
                //     std::string targetType = types[targetIndexLocal + startIdx];
                //     for(int sourceIndexLocal = 0; sourceIndexLocal < sourceWorkload; sourceIndexLocal++){
                //         std::string sourceType = types[sourceIndexLocal + sourceStartIdx];
                //         // if there is at least an interaction between the two types, it means that the virtual outputs are present
                //         if(mappedVirtualOutputsVectors.contains(std::make_pair(sourceType,targetType))){
                //             for(auto virtualOutputPair : mappedVirtualOutputsVectors[std::make_pair(sourceType,targetType)]){
                //                 typeComputations[targetIndexLocal]->setInputVinForType(sourceType, rankVirtualInputsBuffer[sourceRank][virtualInputPosition], virtualOutputPair.first);
                //                 virtualInputPosition++;
                //             }
                //         }
                //     }
                // }
                for (int sourceRank = 0; sourceRank < numProcesses; sourceRank++){
                    int sourceStartIdx = sourceRank * workloadPerProcess;
                    int sourceWorkload;
                    if(sourceRank == (numProcesses-1)){
                        sourceWorkload = types.size() - (sourceRank*workloadPerProcess);
                    } else {
                        sourceWorkload = workloadPerProcess;
                    }
                    int currentVirtualInputIndex = 0;
                    for(int sourceTypeIndex = sourceStartIdx; sourceTypeIndex < sourceStartIdx + sourceWorkload; sourceTypeIndex++){
                        std::string sourceType = types[sourceTypeIndex];
                        for(int targetTypeIndex = startIdx; targetTypeIndex < startIdx + finalWorkload; targetTypeIndex++){
                            std::string targetType = types[targetTypeIndex];
                            std::pair<std::string,std::string> keyTypes = std::make_pair(sourceType, targetType);
                            for(uint localVirtualInputIndex = 0; localVirtualInputIndex < mappedVirtualInputsVectors[keyTypes].size(); localVirtualInputIndex++){
                                std::pair<std::string, std::string> virtualInputPair = mappedVirtualInputsVectors[keyTypes][localVirtualInputIndex];
                                std::pair<std::string, std::string> virtualOutputPair = mappedVirtualOutputsVectors[keyTypes][localVirtualInputIndex];
                                std::string sourceNode = virtualOutputPair.first;
                                std::string targetNode = virtualInputPair.second;
                                std::tuple<std::string, std::string, std::string, std::string> interactionKey = std::make_tuple(sourceNode, targetNode, sourceType, targetType);
                                // typeComputations[targetTypeIndex - startIdx]->setInputVinForType(sourceType, rankVirtualInputsBuffer[sourceRank][virtualInputIndex], virtualInputPair.first);
                                if(interactionBetweenTypesFinerMap[interactionKey].contains(iterationInterType)){
                                    if(sourceType == targetType){
                                        if(sameTypeCommunication) typeComputations[targetTypeIndex - startIdx]->setInputNodeValue(virtualInputPair.first, rankVirtualInputsBuffer[sourceRank][currentVirtualInputIndex]);
                                    } else {
                                        typeComputations[targetTypeIndex - startIdx]->setInputNodeValue(virtualInputPair.first, rankVirtualInputsBuffer[sourceRank][currentVirtualInputIndex]);
                                    }
                                }
                                currentVirtualInputIndex++;
                            }
                        }
                    }
                }
            }
        }

        // delete the virtual inputs vector of arrays
        for(int i = 0; i < numProcesses; i++){
            delete[] rankVirtualInputsBuffer.at(i);
        }
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
