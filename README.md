# c2c-sepia
Cell-to-cell communication with enriched pathways project, using scRNA-seq data.

A framework to work with scRNA-seq profiles is presented here, the algorithm will take as an input:
- The scRNA-seq profiles and log2fold-change
- The type for every cell.

## DEPENDENCIES
- armadillo
- mlpack(latest)
- FUTURE repast HPC
- FUTURE RELOGO

## BUILD
```bash
mkdir build
cmake .
make
```

## USAGE
```bash
./build/c2c-sepia --fMETAPATHWAY [metapathway].tsv --fLogfoldPerCell [logfoldPerCell].tsv --dirCellInteraction [celltypesInteractionFolder]
```

The options are the following:
- **--help**  => print help section
- **--fMETAPATHWAY** (obligatory for now, but it will be a default with the metapathway in the resource folder) => metapathway filename used in the algorithm
- **--fLogfoldPerCell** (obligatory) => logFoldChange per cell matrix
- **--dirCellInteraction** => logfoldPerCell matrix filename
- **--ensembleGeneNames**" (no parameter) => use ensemble gene names, since the metapathway used in resources have entrez_ids, a map will be done from ensemble to entrez, the map is available in resources, if the metapathway is consistent with the data used for the log-fold changes and the interactions (the have the same gene names), this parameter should not be used. Only use in case of external data sources that have ensemble ids and the metapathway used is the one in the resources.
- **--sameCellCommunication** (no parameter) => "use same cell communication, since it is not permitted as the standard definition of the model, this adds a virtual node for the same cell type")
    

## INPUT SCHEMA
**metapathway.tsv**

start \t end \t weight

[gene1] \t [gene2] \t [0.something]

...



**logfoldPerCell.tsv**

\t cell1 \t cell2 \t ... \t cellN 

gene1 \t [lfc_cell1:gene1] \t [lfc_cell2:gene1] \t ... \t [lfc_cellN:gene1]

gene2 \t [lfc_cell1:gene2] \t [lfc_cell2:gene2] \t ... \t [lfc_cellN:gene2]

...


The files in the folder that contains the cell interactions should have the following schema
**celltypesInteraction.tsv**

startCell \t geneLigand \t endCell \t geneReceptor \t weight

[cell1] \t [geneLigand] \t [cell2] \t [genereceptor] \t [0.something]

...



**nsubcelltypes.txt**

cell1

cell3

...



## things to consider
The map file that is used for the mapping of ensemble gene names to entrez ids (required if the metapathway in the __resourse__ folder are used) can b e found in the metapathways folders themselves. The one used in the tool is found at __resources/graphs/metapathwayReactome/nodes.txt__


Under development

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
