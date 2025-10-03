# Comparative single cell database app docker file

# get shiny serveR and a version of R from the rocker project
FROM rocker/shiny:4.0.5

# system libraries
# Try to only install system libraries you actually need
# Package Manager is a good resource to help discover system deps
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libglpk-dev \
    libxml2-dev \
    libbz2-dev \
    liblzma-dev \
    ncbi-blast+ \
    git

# install R packages required 
# Change the packages list to suit your needs

#RUN R -e 'install.packages("ggplot2", repos = "https://packagemanager.rstudio.com/cran/latest")'
#RUN R -e 'install.packages("scales", repos = "https://packagemanager.rstudio.com/cran/latest")'

RUN R -e 'install.packages(c(\
              "ape", \
              "BiocManager", \
              "circlize", \
              "data.table", \
              "devtools", \
              "dplyr", \
              "DT", \
              "eulerr", \
              "ggplot2", \
              "ggpubr", \
              "ggrepel", \
              "gplots", \
              "kableExtra", \
              "Matrix", \
              "patchwork", \
              "RColorBrewer", \
              "R.utils", \
              "remotes", \
              "reshape2", \
              "scales", \
              "shiny", \
              "shinycssloaders", \
              "shinydashboard", \
              "shinyWidgets", \
              "stringr", \
              "tidytree", \
              "vctrs", \
              "viridis", \
              "XML", \
              "zoo" \
            ), \
            #repos="https://packagemanager.rstudio.com/cran/__linux__/focal/latest"\
            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23"\
          )'

RUN R -e 'BiocManager::install(c(\
              "ComplexHeatmap", \
              "InteractiveComplexHeatmap", \
              "treeio" \
          ))'

RUN R -e 'BiocManager::install("ggtree")'

RUN R -e 'install.packages("ggplot2", repos = "https://packagemanager.rstudio.com/cran/latest")'

RUN R -e 'install.packages("ggiraph", repos = "https://packagemanager.rstudio.com/cran/latest")'

RUN R -e 'devtools::install_github( \
              repo="tanaylab/metacell", \ 
              force=TRUE \
          )'


# run app
CMD ["/usr/bin/shiny-server"]
