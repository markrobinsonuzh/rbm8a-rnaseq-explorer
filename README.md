# rbm8a-rnaseq-explorer

In this repository, you will find code and R objects to support the R/Shiny app for exploring rbm8a knockout RNA-seq experiments interactively.

Please see the corresponding preprint:  
Kocere, A.; Chiavacci, E.; Mendez-Acevedo, K. M.; Soneson, C.; Hiltabidle, M. S.; Raghunath, A.; Macgowan, J. S.; Shavit, J. A.; Panakova, D.; Williams, M. L. K.; Robinson, M. D.; Mosimann, C.; Burger, A.: *TAR Syndrome-associated Rbm8a deficiency causes hematopoietic defects and attenuates Wnt/PCP signaling*. bioRxiv doi: [https://doi.org/10.1101/2023.04.12.536513](https://doi.org/10.1101/2023.04.12.536513). Posted April 20, 2023.

Instructions follow below for viewing the dataset locally via R.

# Local Setup

If R is not already installed, find your closest mirror [here](https://cran.r-project.org/mirrors.html) and follow the instructions for installation; if using RStudio, you can download it from [here](https://posit.co/download/rstudio-desktop/))

Download the repository, either from the command line with a command like:

```
git clone https://github.com/markrobinsonuzh/rbm8a-rnaseq-explorer.git
```

or, download a snapshot of the repository using the green `Code` button above (e.g., click on `Download ZIP` and then unzip the files of the repository in a reasonable location locally).

Then, start an R session and set the current working directory to where `rbm8a-rnaseq-explorer` was placed. (TIP: if using RStudio, you can double-click on the `rbm8a-rnaseq-explorer.Rproj` file and an R session will be started with the directory correctly set.)

If not already installed, you will need to install the following packages (cut-and-paste into the console):

```
install.packages("BiocManager", quiet = TRUE)
pkgs <- c("shiny","ggplot2","ggrepel",
          "GGally", "dplyr", "tidyr",
          "Gviz", "ComplexHeatmap", "reshape2",
          "circlize", "gridExtra", "rtracklayer")
BiocManager::install(pkgs)
```

To run the app, run:

```
source("app.R")
```

Note that for the coverage tracks to be displayed, you will need to download 7 BW (bigwig) files from [here](https://doi.org/10.5281/zenodo.7680271) and place them in the `starbigwig` directory.


# About (the panels of the app)

See [about_app.md](about_app.md).
