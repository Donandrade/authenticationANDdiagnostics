# authenticationANDdiagnostics
Analysis for the selection of genetic targets aiming at product authentication and diagnosis in a broader context.

HRM Analysis

This repository contains an R script for performing High-Resolution Melting (HRM) analysis on raw data. The script processes HRM data, calculates dissimilarity matrices, performs hierarchical clustering, and generates various plots.

## Dependencies

The script requires the following R packages. You can install them using the commands below:

```r
install.packages("tidyr")
if (!require("processx")) install.packages("processx")
install.packages("MBmca")

# Additionally, load the necessary libraries
library(ggplot2)
library(factoextra)
library(dendextend)
library(cluster)
library(ape)
library(circlize)
library(qpcR)
library(tidyverse)
library(plotly)
library(MBmca)
library(RColorBrewer)
```

Script Overview
1. Data Loading and Filtering
Reads the HRM raw data from a file located at ../data/M144_Raw_Data.txt.
Filters the data based on specific temperature ranges (77°C to 84°C) and adjusts column names for consistency.
2. Melt Curve Analysis
Performs melt curve analysis to normalize the data and extract fluorescence values. This analysis is used to adjust for baseline shifts and to compute fluorescence intensities.
3. Dissimilarity Calculation
Computes a dissimilarity matrix considering genotype confidence percentage. This matrix quantifies the differences between samples based on their HRM profiles.
Saves the dissimilarity matrix to dissimilaridade.txt.
4. Hierarchical Clustering
Performs hierarchical clustering on the dissimilarity matrix using average linkage method.
Generates a dendrogram plot to visualize the clustering results.
Saves the dendrogram plot as M144.pdf.
5. Cluster Identification and Visualization
Identifies clusters in the hierarchical clustering results and assigns colors to each cluster.
Generates a line plot of the HRM data, color-coded by cluster, to show the fluorescence data across different temperature ranges.
Saves the line plot as df_M144.pdf.
6. Cluster Plot
Plots a colored dendrogram with labels indicating different clusters. The colors represent different clusters and provide a visual indication of how samples are grouped based on their dissimilarities.
Plot Details
Dendrogram Plot (M144.pdf):

Shows the hierarchical clustering results as a dendrogram.
The plot illustrates how samples are grouped into clusters based on their dissimilarity scores.
Line Plot (df_M144.pdf):

Displays the fluorescence data against temperature for each sample.
Lines are color-coded according to the clusters identified in the hierarchical clustering.
Usage
Place your HRM raw data file in the ../data/ directory with the filename M144_Raw_Data.txt.

Run the R script to process the data, calculate dissimilarity, perform clustering, and generate plots.

Check the output files:

dissimilaridade.txt for the dissimilarity matrix which includes the effect of genotype confidence percentage.
M144.pdf for the hierarchical clustering dendrogram.
df_M144.pdf for the line plot of HRM data.
Notes
Ensure that all necessary R packages are installed.
The script assumes that the HRM data file is formatted correctly and located in the specified directory.
Adjust the script as needed for different data formats or analysis requirements.
