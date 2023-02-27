# rbm8a-rnaseq-explorer

In this repository, you will find code and R objects to support the R/Shiny app for exploring rbm8a knockout RNA-seq experiments interactively.

TODO: citation (Kocere et al.)

The app is available from [here](http://imlspenticton.uzh.ch:3838/mosimann_p2452/). 

Instructions follow below for viewing the dataset locally via R.

# Local Setup

First, download the repository, either from the command line with a command like:

```
git clone https://github.com/markrobinsonuzh/rbm8a-rnaseq-explorer.git
```

or, download a snapshot of the repository using the green `Code` button above (e.g., click on `Download ZIP` and then unzip the files of the repository in a reasonably location locally).

Then, start an R session and set the current working directory to whereever `rbm8a-rnaseq-explorer` was placed. (TIP: you can double-click on the `rbm8a-rnaseq-explorer.Rproj` and an R session will be started with the directory correctly set.)

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

# About (the panels of the app)

## PCA
The PCA is based on the rlog-transformed values from DESeq2. You can select which pair of principal components to show. The axes also indicate the fraction of the variance explained by each principal component. 

## CPM plot
The CPMs shown in this plot are calculated by edgeR. You can search for a gene of interest in the `Selected gene` input box. You can provide either Ensembl IDs (without the final version number) or gene symbols. If the search matches multiple genes, they will all be shown. The search is case-insensitive. Hovering over the points shows the corresponding CPM values. 

## Coverage plot
The coverage plot is obtained from the BAM files generated by STAR. The gene model is shown in blue under the coverage tracks. There are two additional gene model tracks (although they are empty for many genes): in green are all features from other genes that fall in the considered genomic region, and in orange are all the transcripts from the gene of interest that were excluded from the reference since they were not expressed highly enough in any of the samples. You can search for a gene of interest in the `Selected gene` input box. You can provide either Ensembl IDs (without the final version number) or gene symbols. If the search matches multiple genes, only one will be shown (thus, it is safer to search by Ensembl ID here, or to make sure via the CPM plot that there is only one gene with a given symbol). The search is case-insensitive. The sample ID and the group are indicated in the left title panel of each track. Note that the coverage range differs between tracks, and that the library size is not taken into account, so the absolute coverage numbers are not necessarily comparable across samples. 

## Volcano plots
The volcano plots correspond to the differential gene expression results obtained by `edgeR`, for the four different contrasts. You can search for a gene of interest in the `Selected gene` input box, which will highlight it in the volcano plots. You can provide either Ensembl IDs (without the final version number) or gene symbols. If the search matches multiple genes, they will all be shown. The search is case-insensitive. 

## DGE result table
The result table corresponds to the differential gene expression results obtained by `edgeR`, for the four different contrasts. The result table can be searched using the search box in the top right, and can be ordered by each column. Note that this search is independent from the `Selected gene` input box, and that in the result table, the search will return also partial matches. 

## DEXSeq result table
The result table corresponds to the differential exon/intron usage results obtained by `DEXSeq`, for the four different contrasts. The result table can be searched using the search box in the top right, and can be ordered by each column. Note that this search is independent from the `Selected gene` input box, and that in the result table, the search will return also partial matches. The `DEXSeq` analysis searches for features (exon bins/introns) that are differentially *used* (included/excluded) between conditions, and thus disregards any overall differences in gene expression. Exons and introns are numbered "from left to right" in terms of genome coordinates, regardless of the strand. 

## SUPPA result table
The result table corresponds to the event-based differential splicing results obtained by `SUPPA`, for the four different contrasts. The result table can be searched using the search box in the top right, and can be ordered by each column. Note that this search is independent from the `Selected gene` input box, and that in the result table, the search will return also partial matches. The different event types are:

- A3: alternative 3' splice site
- A5: alternative 5' splice site
- AF: alternative first exon
- AL: alternative last exon
- MX: mutually exclusive exons
- RI: retained intron
- SE: skipped exon

## DEXSeq DTU result table
The result table corresponds to the differential transcript usage (DTU) results obtained by applying DEXSeq to transcript counts (rather than exon/intron bin counts). Note that all genes that only have one (retained) isoforms will have NA results in this analysis, since they can not show differential usage of their isoforms. 

## Troubleshooting
If the coverage plot tab gives an error like: 
```
> Invalid chromosome identifier 'X'
Please consider setting options(ucscChromosomeNames=FALSE) to allow for arbitrary chromosome identifiers.
```

then delete the selected gene in the search box and type it in again.  
