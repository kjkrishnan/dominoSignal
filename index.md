---
output: github_document
---

<!-- index.md is generated from index.Rmd. Please edit that file -->



# Introducing domino2

domino2 is an updated version of the original [domino](https://github.com/Elisseeff-Lab/domino) R package published in Nature Biomedical Engineering in [Computational reconstruction of the signalling networks surrounding implanted biomaterials from single-cell transcriptomics](https://doi.org/10.1038/s41551-021-00770-5). domino2 is a tool for analysis of intra- and intercellular signaling in single cell RNA sequencing data based on transcription factor activation and receptor and ligand linkages.

## Installation

domino2 is undergoing active development where aspects of how data is used, analyzed, and interpreted is subject to change as new features and fixes are implemented. v0.2.1 of domino2 serves as the first stable development version during these active updates for reproducible usage (see our [News](NEWS.md) page for more information on changes).

This version is currently hosted on the [FertigLab GitHub](https://github.com/FertigLab) as [branch v0.2.1](https://github.com/FertigLab/domino_development/tree/v0.2.1) of the [domino_development repository](https://github.com/FertigLab/domino_development) forked from the primary repository hosted on the [Elisseeff-Lab GitHub](https://github.com/Elisseeff-Lab/domino), and can be installed using the remotes package.


```r
if(!require(remotes)){
    install.packages('remotes')
}
remotes::install_github('FertigLab/domino_development@v0.2.1')
```

## Usage Overview

Here is an overview of how domino2 might be used in analysis of a single cell RNA sequencing data set:

1. Transcription factor activation scores are calculated (we recommend using [pySCENIC](https://pyscenic.readthedocs.io/en/latest/), but other methods can be used as well)
2. A ligand-receptor database is used to map linkages between ligands and receptors (we recommend using [CellphoneDB](https://www.cellphonedb.org/), but other methods can be used as well).
3. A domino object is created using counts, z-scored counts, clustering information, and the data from steps 1 and 2.
4. Parameters such as the maximum number of transcription factors and receptors or the minimum correlation threshold (among others) are used to make a cell communication network
5. Communication networks can be extracted from within the domino object or visualized using a variety of plotting functions

Please see the [Getting Started](vignette("domino2")) page for an example analysis that includes all of these steps in detail, from downloading and running [pySCENIC](https://pyscenic.readthedocs.io/en/latest/) to building and visualizing domino results. Other articles include further details on [plotting functions](vignette("plotting_vignette")) and the structure of the [domino object](vignette("domino_object_vignette")).