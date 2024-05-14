# R-utils

A collection of utility scripts in R for data analysis and visualization to serve as a reference for lab members.

## How to use

In general, the scripts are prepared to be downloaded and manually edited. There is usually an `INPUTS` section at the beginning of the script where you can specify the input files and parameters. The scripts are prepared to be run from RStudio, but you can also run them in the terminal with `Rscript <scriptname>`.

## Alternative splicing plots

- `AS_boxplot.R`: Boxplot of PSI values for alternative splicing events. Includes a Wilcoxon test for differential splicing between all of the groups.
- `AS_kmeans.R`: K-means clustering of PSI values for alternative splicing events.
- `AS_PCA-heatmap-UMAP.R`: PCA, heatmap, and UMAP plots for alternative splicing data. Allows to filter events by overall variance. Specific event and sample-lists can be provided.


## Gene expression plots

- `GE_CGGA.R`: PCA, heatmap and UMAP plots for gene expression data. Used for the CGGA dataset. 
- `GE_plots.R`: general script with boxplots with tests, PCA, heatmap and UMAP plots. Allows to subset specific genes and samples, as well as to filter genes by overall variance.