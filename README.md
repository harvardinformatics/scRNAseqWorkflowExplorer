# Comparing outputs from scRNA-seq workflows
This repostiory contains a Shiny app the purpose of which is to generate visualizations for comparing the outputs of multiple scRNA-seq workflows run on a single scRNA-seq sample library. The app has been designed to take as input the outputs from two Snakemake workflows: 

* A workflow that generates scRNA-seq Seurat objects stored in R *rds* files, including the endpoints of 6 pre-processing tool combinations, and 11 additional rds files for intermediate steps in the workflow tree; this workflow can be found [here](https://github.com/harvardinformatics/scRNAseq-preprocessing).
* A workflow that performs downsampling and reclustering on Seurat objects, compares clustering from the downsampled data to the original full data to quantify cluster stability,and writes tabular data summarizing stability; the workflow can be found [here](https://github.com/harvardinformatics/scRNAseqClusterDownsample).

**NOTE:** The first of these workflows generates annotations in the object metadata that are used by this app, so using this app with data generated outside of the above Snakemake workflows will likely require generating those annotations.

**DO NOT** use this app to compare files for two different sequencing libraries as such comparisons are meaningless.


## Running the app
For lightweight analyses for a subset of the files produced by, running the app locally on a laptop should work, so long as the number of cells in the Seurat objects are not too large. Attempting to load all 17 Seurat objects produced by the Snakemake workflow requires > 18Gb of memory, i.e. the typical maximum available on a Macbook pro. A more robust solution is to launch the app on an HPC cluster, such as where you ran the Snakemake workflows. This requires having an option for launching Rstudio Server on a remote machine in a way that a browser window will open up the app locally. For users of Harvard's Cannon cluster, the easiest way to do this is to launch the server via [Open OnDemand](https://docs.rc.fas.harvard.edu/kb/virtual-desktop/). Once you've launched the OOD job, and open Rstudio, follow the steps below. 

### R packages
Wherever you run the app, you will need to make sure that the relevant R packages are installed. in the R console, you will need to install them like this:


```bash
install.packages(
  c("Seurat", "SeuratObject", "shiny", "tidyverse",
    "ComplexUpset", "ggrepel", "ape", "Matrix"),
  lib = .libPaths()[1],
  repos = "https://cloud.r-project.org"
)
```

### Data files
Put your *.rds* files and downsampling (.tsv)  files output by the Snakemake workflow in the *data/* directory. Next,



### Launch the app
There are two main ways to launch the app. One can launch it from a terminal window.

From inside the `shiny/` directory:

```bash
Rscript -e "shiny::runApp('.')"
```

Or from anywhere using the full path:

```bash
Rscript -e "shiny::runApp('/PATH/TO/shiny')"
```

Alternatively, you can launch it from the R (typically Rstudio) console

```bash
setwd("/path/to/app")
shiny::runApp(".")
```

