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
./build/c2c-sepia <???>
```

Under development

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
