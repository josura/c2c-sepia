#include <iostream>
#include "Computation.h"
#include "WeightedEdgeGraph.h"
#include "utilities.h"

void printHelp(){
    std::cout << "usage: ./c2c-sepia <metapathway>.tsv <logfoldPerCell>.tsv [<subcelltypes>.txt] [<celltypesInteraction>.tsv]\nFILE STRUCTURE SCHEMA:\nmetapathway.tsv\nstart\tend\tweight\n<gene1>\t<gene2>\t <0.something>\n...\n\n\nlogfoldPerCell.tsv\n\tcell1\tcell2\t...\tcellN\ngene1\t<lfc_cell1:gene1>\t<lfc_cell2:gene1>\t...\t<lfc_cellN:gene1>\ngene1\t<lfc_cell1:gene2>\t<lfc_cell2:gene2>\t...\t<lfc_cellN:gene2>\n...\n\n\ncelltypesInteraction.tsv\nstartCell:geneLigand\tendCell:geneReceptor\tweight\n<cell1:geneLigand>\t<cell2:genereceptor>\t <0.something>\n...\n\n\nsubcelltypes.txt\ncell1\ncell3\n..."<<std::endl;
    std::cout << "LEGEND:\n <> := placeholder for the name of the file\n[] := optional\n{} := at least one"<<std::endl;
}

int main(int argc, char** argv ) {
    WeightedEdgeGraph* g1 = new WeightedEdgeGraph(4);
    WeightedEdgeGraph* g2 = new WeightedEdgeGraph();
    std::cout << " testing "<< g2->getNumEdges()<<std::endl;
    std::cout << *g2 ;
    g1->addEdge(1,2,0.1)->addEdge(1,3,0.4);
    *g2 = *g1;
    WeightedEdgeGraph graphtest = *g1;
    std::cout << " testing "<< g1->getNumEdges()<<std::endl;
    std::cout << *g1 ;
    std::cout << " testing "<< g2->getNumEdges()<<std::endl;
    std::cout << *g2 ;

    std::vector<std::string> nodeNames{"node1","node2","node3","node4","node5"};
    std::vector<double> nodeValues{0.3,4.1,3.8,8.2,9.5};
    WeightedEdgeGraph g4 = WeightedEdgeGraph(nodeNames);
    WeightedEdgeGraph g5 = WeightedEdgeGraph(nodeNames,nodeValues);
    WeightedEdgeGraph g6;
    g6.assign(g5);

    g5.addEdge("node1","node3",2.3)->addEdge("node3","node2",0.3);

    std::cout << " testing "<< g4.getNumEdges()<<std::endl;
    std::cout << g4 ;
    std::cout << " testing "<< g5.getNumEdges()<<std::endl;
    std::cout << g5 ;
    //std::cout << " testing "<< g6.getNumEdges()<<std::endl;
    //std::cout << g6 ;

    std::string thisCellType = "testCell";
    std::vector<double> input{0.1,0.3,0.5,0.62,0.34,0.87}; 
    //Matrix<double> _W=Matrix<double>::createRandom(6, 6); 
    std::vector<double> matrixVector{0,0.2,0.4,0.5,0.23,0.67,
                                        0.2,0,0.4,0.5,0.23,0.67,
                                        0.3,0.2,0,0.5,0.23,0.67,
                                        0.43,0.2,0.4,0,0.23,0.67,
                                        0.54,0.2,0.4,0.5,0,0.67,
                                        0.25,0.2,0.4,0.5,0.23,0,};
    Matrix<double> _W = Matrix<double>(matrixVector,6,6);
    std::vector<std::string> metapathwayNames{"testGene1","testGene2","testGene3","testGene4","testGene5","testGene6"};

    const std::vector<std::string> cellTypes{"testCell","testCell2","testCell3","testCell4"};
    const std::vector<std::pair<std::string, std::string>> virtualInputEdges{{"v-in:testCell2","testGene2"},
                                                                                {"v-in:testCell4","testGene2"},
                                                                                {"v-in:testCell4","testGene3"},
                                                                                {"v-in:testCell4","testGene6"}
                                                                                };
    std::vector<double> virtualInputEdgesValues = {0.4,0.5,0.7,0.2};
    std::vector<std::pair<std::string, std::string>> virtualOutputEdges{{"testGene2","v-out:testCell2"},
                                                                                        {"testGene4","v-out:testCell3"}
                                                                                    };
    std::vector<double> virtualOutputEdgesValues = {0.4,0.4};

    
    auto c1  = new Computation(thisCellType,input,_W,metapathwayNames);
    c1->augmentMetapathway(cellTypes);
    c1->addEdges(virtualInputEdges,virtualInputEdgesValues);
    c1->addEdges(virtualOutputEdges,virtualOutputEdgesValues);
    auto perturbation = c1->computeAugmentedPerturbation();
    c1->updateInput(perturbation,true);
    auto computationInputAfter = c1->getInputAugmented();

    std::string filename,celltypesFilename,celltypesInteractionFilename,cellLogFoldMatrixFilename;
    if(argc < 2){
        printHelp();
        return 0;
    }else{
        if(argc >= 2)
            filename = argv[1];
        if(argc >= 3){
            cellLogFoldMatrixFilename = argv[2];  // first row is cell types, so default means all the celltypes are used
        }
        if(argc >= 4){
            celltypesInteractionFilename = argv[3];
        }
        if(argc >= 5){
            celltypesFilename = argv[4]; //all the celltypes or a subset will be used
        }
    }
    auto namesAndEdges = edgesFileToEdgesListAndNodesByName(filename);
    
    return 0;
}