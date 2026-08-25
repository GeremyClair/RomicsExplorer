# RomicsExplorer
Shiny app designed to explore bulk omics datasets in the Romics_object format. The repository contains two versions of the app: one that allows loading any romics_object self-contained dataset, and another that is preloaded with an example Bacillus cereus dataset (see RomicsProcessor preprint for details).

# RomicsProcessor
RomicsProcessor is an R package for analyzing bulk, single-cell, and spatial omics datasets.
The package provides a structured R object (`romics_object`) to store data, metadata, and complete processing history,
enabling reproducible and FAIR-compatible data analysis. RomicsProcessor enables complete analytical traceability through
UUID-based object tracking and automatic step logging, supporting reusable analytical pipelines and reliable research reproducibility.

### Prerequisites
Note: it is possible that packages not listed below are required; please read the R error message carefully to identify those.

```R
install.packages(“devtools”)
install.packages("shinyapp")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("ComplexHeatmap")
```

### Install RomicsProcessorfrom GitHub
```R
devtools::install_github("PNNL-Comp-Mass-Spec/RomicsProcessor")
```

### Start the shiny app
After decompressing all files in a single folder go into RStudio and open the shiny app code file. Then press on the green play button to start the app.

## System Requirements

- **R**: ≥ 4.5.1 (tested with R 4.5.1, recommended for case study reproducibility)
- **Memory**: For large datasets (>100k samples), recommend ≥32GB RAM
 
## Getting Started - Case Study Examples

RomicsProcessor includes a comprehensive case study demonstrating the complete workflow:

## Cite the code

To cite the package, please use the following DOI:
[![citation(https://www.biorxiv.org/content/10.64898/2026.07.09.737600v1)

## Authors

**Geremy Clair** (main contact) - geremy.clair@pnnl.gov \
Harsh Bhotika - harsh.bhotika@pnnl.gov \
Brittney Gorman - brittney.gorman@pnnl.gov

Written for the Department of Energy (PNNL, Richland, WA) \
Website: https://omics.pnl.gov/ or https://panomics.pnnl.gov/ \
E-mail: proteomics@pnnl.gov

## License

RomicsProcessor is licensed under the 3-Clause BSD License; 
you may not use this file except in compliance with the License.  You may obtain 
a copy of the License at https://opensource.org/licenses/BSD-3-Clause

Copyright 2026 Battelle Memorial Institute

