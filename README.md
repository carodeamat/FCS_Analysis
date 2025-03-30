# FCS Analysis
Analysis of high-parameter flow cytometry data for the identification and comparison of cell populations through unsupervised clustering.

## Input data
The workflow requires FCS files from population of interest. \
That is, after excluding doublets, dead cells, and any other population using FlowJo. \
The population is exported from FlowJo with compensated parameters.

This workflow include three main steps:
1. Data transformation.
2. Clustering using SOM or Phenograph.
3. Annotation and visualization of results.


## Data transformation
The `01_TransformData.Rmd` file contains the workflow code to read the FCS files and a metadata file. \
The metadata file is a csv file which contains in the first column the file name of the FCS files that will be analyzed. \
The data will be transformed using the arcsinh algorithm with adjustable cofactors. \
Once the data is transformed, two outputs will be saved:
- Expression matrix with transformed values.
- FCS files with transformed values.

## Clustering
There are two main methods for clustering flow cytometry data, FlowSOM and Phenograph. \
Each of them has their own strengths and limitations. \
`02a_SOM_clustering.Rmd` is the workflow to perform Flowsom clustering. \
`02b_Pheno_clustering.Rmd` is the workflow to perform Phenograph clustering. \
Both require the expression matrix generated after transformation and the metadata file. \
tSNE and UMAP plots are generated to evaluate clustering parameter. \
The output is an rds file containing a list with expresion matrix and clustering information. \

## Annotation and visualization
The `03_Annotation_Visualization.Rmd` file has the workflow to evaluate the expression of markers in each cluster as heatmaps and 2D tSNE or UMAP plots. \
The clusters can be renamed after marker evaluation. \
Descriptive bar plots and density plots can also be genreated. \
The output is and rds file with the annotated data and svg plots.
