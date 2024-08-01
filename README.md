# authenticationANDdiagnostics
Analysis for the selection of genetic targets aiming at product authentication and diagnosis in a broader context.

HRM Analysis

This repository contains an R script (`hrm_analysis/m144.R`) for performing High-Resolution Melting (HRM) analysis on raw data. The script processes HRM data, calculates dissimilarity matrices, based on Genotype Confidence Percentage (**GCP**) performs hierarchical clustering, and generates various plots.

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

## Script Overview
### Data Loading and Filtering
Reads the HRM raw data from a file located at `../data/M144_Raw_Data.txt`.
Filters the data based on specific temperature ranges (**77°C** to **84°C**) and adjusts column names for consistency.
### Melt Curve Analysis
Performs melt curve analysis to normalize the data and extract fluorescence values. This analysis is used to adjust for baseline shifts and to compute fluorescence intensities.
### Dissimilarity Calculation
- Computes a dissimilarity matrix considering genotype confidence percentage. This matrix quantifies the differences between - - samples based on their HRM profiles.



### Calculation of Genotype Confidence Percentage (GCP)

1. **Calculation of ![Formula \( S_{rt} \)](https://latex.codecogs.com/svg.latex?S_{rt}):**

   The formula is:
   ![Formula \( S_{rt} \)](https://latex.codecogs.com/svg.latex?S_{rt}%20=%201.05^{-0.02%20\sum_{i=a}^{z}\left(f_{ri}%20-%20f_{ti}\right)^{2}})

   **Explanation:**
   - ![Formula \( S_{rt} \)](https://latex.codecogs.com/svg.latex?S_{rt}) is a measure of similarity.
   - The base of the exponentiation is
     ![Base](https://latex.codecogs.com/svg.latex?1.05)
   - The coefficient that adjusts the influence of the sum of squared differences is
     ![Coefficient](https://latex.codecogs.com/svg.latex?-0.02)
   - The summation spans an index range ![Formula](https://latex.codecogs.com/svg.latex?{i}) from ![Formula](https://latex.codecogs.com/svg.latex?{a}) to ![Formula](https://latex.codecogs.com/svg.latex?{z}): 
     ![Summation](https://latex.codecogs.com/svg.latex?\sum_{i=a}^{z})
   - The term ![Term](https://latex.codecogs.com/svg.latex?\left(f_{ri}%20-%20f_{ti}\right)^{2}) is the square of the difference between the fluorescence values ![Term](https://latex.codecogs.com/svg.latex?f_{ti}%20) and ![Term](https://latex.codecogs.com/svg.latex?f_{ri}%20) for each index ![Formula](https://latex.codecogs.com/svg.latex?{i})
   - The result of the summation is multiplied by \( -0.02 \) and used as the exponent for \( 1.05 \):
     ![Summation Result](https://latex.codecogs.com/svg.latex?1.05^{-0.02%20\sum_{i=a}^{z}\left(f_{ri}%20-%20f_{ti}\right)^{2}})

2. **Calculation of![Formula \( S_{rt} \)](https://latex.codecogs.com/svg.latex?D_{rt}):**

   The formula is: ![Formula \( D_{rt} \)](https://latex.codecogs.com/svg.latex?D_{rt}%20=%20(1%20-%20S_{rt})%20*%20100)

   **Explanation:**
   - ![Formula \( S_{rt} \)](https://latex.codecogs.com/svg.latex?D_{rt}) is the measure of dissimilarity.
   - It is calculated as the difference between 1 and ![Formula \( S_{rt} \)](https://latex.codecogs.com/svg.latex?S_{rt}), multiplied by 100 to convert it into a percentage:
     ![Calculation of \( D_{rt} \)](https://latex.codecogs.com/svg.latex?D_{rt}%20=%20(1%20-%20S_{rt})%20*%20100)
   - The higher the value of ![Formula \( S_{rt} \)](https://latex.codecogs.com/svg.latex?D_{rt}), the greater the dissimilarity between the fluorescence profiles.

The ![Formula \( S_{rt} \)](https://latex.codecogs.com/svg.latex?D) matrix is saved to `dissimilaridade.txt` file.

### Hierarchical Clustering
- Performs hierarchical clustering on the dissimilarity matrix using average linkage method.
- Generates a dendrogram plot to visualize the clustering results.
- Saves the dendrogram plot as `M144.pdf`.
### Cluster Identification and Visualization
- Identifies clusters in the hierarchical clustering results and assigns colors to each cluster.
- Generates a line plot of the HRM data, color-coded by cluster, to show the fluorescence data across different temperature ranges.
- Saves the line plot as `df_M144.pdf`.
### Cluster Plot
Plots a colored dendrogram with labels indicating different clusters. The colors represent different clusters and provide a visual indication of how samples are grouped based on their dissimilarities.

## USAGE
Place your HRM data file in the `../data/ directory`. The example file used is named M144_Raw_Data.txt. You can use the same name for your file, or if you choose a different name, make sure to update the filename on line 22 of the code (`hrm_analysis/m144.R`).

Run the R script to process the data, calculate dissimilarity, perform clustering, and generate plots.

## INPUT
An example of the input file is in the data directory. This dataset corresponds to the first derivative of fluorescence obtained from a real-time PCR run on the LightCycler® equipment.

## OUTPUT

- `dissimilaridade.txt` for the dissimilarity matrix which includes the effect of genotype confidence percentage.
- `M144.pdf`:
Shows the hierarchical clustering results as a dendrogram. The plot illustrates how samples are grouped into clusters based on their dissimilarity scores. For now, we suggest you check this file to select the number of K-means clusters and then choose the number of clusters. In our example, we selected `k=2` (see the `k` on line 122 of the script `hrm_analysis/m144.R`).
- `df_M144.pdf`:  for the line plot of the first derivative plot (`-dF/dT`):
Displays the `-dF/dT` of fluorescence data against the shift temperature for each sample. This step will require future optimization, as it currently depends on user adjustments. Therefore, be attentive and choose the best temperature shift for your data.

## NOTES
- Ensure that all necessary R packages are installed.
- The script assumes that the HRM data file is formatted correctly and located in the specified directory.
- Adjust the script as needed for different data formats or analysis requirements.
