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
./build/c2c-sepia <metapathway>.tsv <logfoldPerCell>.tsv [<subcelltypes>.txt] [<celltypesInteraction>.tsv]  
```

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



**celltypesInteraction.tsv**

startCell \t geneLigand \t endCell \t geneReceptor \t weight

[cell1] \t [geneLigand] \t [cell2] \t [genereceptor] \t [0.something]

...



**nsubcelltypes.txt**

cell1

cell3

...




Under development

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
