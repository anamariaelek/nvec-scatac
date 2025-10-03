# metacell-downstream-functions

Set of functions for downstream analysis and visualization of scRNA data processed with Metacell.

Contents:

* `Cross_species_functions.R`: functions for cross-species analysis
* `Downstream_functions.R`: analysis (`sca`) and plotting (`scp`) functions.
* `Export_functions.R`: ?
* `Modified_functions.R`: custom versions of functions from the `metacell` package
* `Tree_functions.R`: functions to build cell type trees, within a single species
* `Gene_module_functions.R`: functions to compute gene modules, based on WGCNA.
* `utils.R`: ?

## Tutorial

There's a simple tutorial [here](https://github.com/sebepedroslab/metacell-examples).

## Install the `metacell` environment

Set up the basic environment, based on `R` 4.0, and include key packages:

```bash
# conda create -n scan
# conda install -c conda-forge r-base=4.0.3 cairo=1.16.0 r-cairo=1.5
conda env create -n scan --file environment.yaml
```

Additional libraries — run within `R`:

```R
install.packages("BiocManager")
install.packages("jsonlite")
BiocManager::install("remotes")
BiocManager::install("tanaylab/metacell")
BiocManager::install("ComplexHeatmap")
# other common packages should be pulled from CRAN/bioconductor (zoo, scales, RColorBrewer, ggplot2...)

# all other packages are available at CRAN
install.packages("tglkmeans", repos=c(getOption("repos"), "https://tanaylab.github.io/repo"))
install.packages("rasterpdf")
install.packages("treemap")
```
