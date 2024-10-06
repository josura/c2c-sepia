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
#include "CustomFunctions.h"
#include "Logger.hxx"
#include "Checkpoint.hxx"



int main(int argc, char** argv) {    
    //program options
    bool sameTypeCommunication=false;
    bool saturation=false;
    bool customSaturation=false;
    bool conservateInitialNorm=false;
    bool undirected = false;
    bool undirectedTypeEdges = false;
    bool resetVirtualOutputs = false;
    bool resumeCheckpoint = false;
    std::string logMode="";
    std::string quantizationMethod = "single";
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
        ("nodeDescriptionFile", po::value<std::string>(), "(string) node description file, used to generate the output description in the case of common graph between types, if not specified no names are used. For an example see in data resources/graphs/metapathwayNew/nodes.tsv")
        ("nodeDescriptionFolder", po::value<std::string>(), "(string) nodes folder, where the files containing the description/nodes for all the graphs are contained, used to read the graph nodes, if not specified the graphs will be built with the edges files(could not contain some isolated nodes) for an example see the folder structure in data/testdata/testHeterogeneousTemporalGraph/nodesDescriptionDifferentStructure")
        ("sameTypeCommunication",po::bool_switch(&sameTypeCommunication),"() use same type communication, since it is not permitted as the standard definition of the model, this adds a virtual node for the same type type")
        ("outputFolder",po::value<std::string>(),"(string) output folder for output of the algorithm at each iteration")
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
        ("saturationTerm",po::value<double>(),"defines the limits of the saturation [-saturationTerm,saturationTerm], default to 1, if saturation is not set, this option is not used, if specified the program will stop the execution")
        ("customSaturationFunction",po::bool_switch(&customSaturation),"use custom saturation function defined in src/CustomFunctions.cxx, if this option is not set, the saturation function will be the default one")
        ("conservateInitialNorm",po::bool_switch(&conservateInitialNorm), "conservate the initial euclidean norm of the perturbation values, that is ||Pn|| <= ||Initial||, default to false")
        ("undirectedEdges",po::bool_switch(&undirected), "edges in the graphs are undirected")
        ("undirectedTypeEdges",po::bool_switch(&undirectedTypeEdges), "edges between types are undirected")
        ("resetVirtualOutputs",po::bool_switch(&resetVirtualOutputs), "reset the virtual outputs to 0 at each iteration, default to false")
        ("virtualNodesGranularity", po::value<std::string>(), "(string) granularity of the virtual nodes, available options are: 'type', 'node'(unstable), 'typeAndNode', default to type")
        ("virtualNodesGranularityParameters", po::value<std::vector<std::string>>()->multitoken(), "(vector<string>) parameters for the virtual nodes granularity, NOT USED for now")
        ("quantizationMethod",po::value<std::string>(), "(string) define the quantization method used to quantize the contact times for the edges between different types, available options are: 'single' and 'multiple'") // aggiungere documentazione
        ("loggingOptions",po::value<std::string>(&logMode),"(string) logging options, available options are: 'all','none', default to all")
        ("savePerformance",po::value<std::string>(&performanceFilename), "(string) output performance (running time, number of total nodes, number of communities, number of total edges) to the defined file, if nothing is specified the performance are not saved")
        ("resumeCheckpoint",po::bool_switch(&resumeCheckpoint), "resume the computation from the last checkpoint, if the checkpoint is not found, the computation will start from the beginning")
    ;

    

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
