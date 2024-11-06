# Script for HRM Data Analysis
# This script processes High-Resolution Melting (HRM) data, calculates dissimilarity between samples,
# performs hierarchical clustering, generates a phylogenetic tree, and visualizes clusters using Plotly.
# Prerequisites: Install packages ("tidyr", "processx", "MBmca") if not already installed.

# install.packages("tidyr")
# if (!require("processx")) install.packages("processx")
# install.packages("MBmca")

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

# Load and Filter Data
data <- read.csv("../hrm_data/M144_Raw_Data.txt", header = TRUE, sep = "\t", dec = ".", check.names = FALSE)
data2 <- c()

# Filter data within temperature range 77°C to 84°C
for(i in 1:nrow(data)) {
  if(data[i,1] > 77 & data[i,1] < 84) {
    tmp <- c(data[i,])
    data2 <- as.data.frame(rbind(as.data.frame(data2), as.data.frame(tmp)))
  }
}

# Adjust column names by removing prefix 'X'
names(data2)[startsWith(names(data2), "X")] <- substring(names(data2)[startsWith(names(data2), "X")], 2)
data <- data2

# Normalize fluorescence data using meltcurve
x <- meltcurve(data %>% select(1, 2), norm = TRUE, span.smooth = 0.05, is.deriv = TRUE, peaklines = TRUE, calc.Area = TRUE)
c <- x[[1]][1]
c <- cbind(c, x[[1]][3])

# Create matrix for normalized fluorescence and temperature data
d <- x[[1]][1]
d <- cbind(d, x[[1]][["Fluo"]])

for(i in 3:ncol(data)) {
  x <- meltcurve(data %>% select(1, i), norm = TRUE, span.smooth = 0.05, is.deriv = TRUE, peaklines = TRUE, calc.Area = TRUE)
  c <- cbind(c, x[[1]][3])
  d <- cbind(d, x[[1]][["Fluo"]])
}

colnames(c) <- colnames(data)
colnames(d) <- colnames(data)

# Filter temperature range for normalized fluorescence
filter <- c()
for(i in 1:nrow(c)) {
  if(c[i,1] > 77 & c[i,1] < 84) {
    tmp <- c(c[i,])
    filter <- as.data.frame(rbind(as.data.frame(filter), as.data.frame(tmp)))
  }
}

colnames(d) <- colnames(data)
normalized <- d[,2:ncol(d)]

# Calculate dissimilarity between samples
matrix <- c()
tmpDiss <- c()
dissimilarity <- c()

for(a in 1:ncol(normalized)) {
  for(j in 1:ncol(normalized)) {
    tmp <- as.matrix((normalized[,a]-normalized[,j])^2)
    matrix <- as.matrix(cbind(matrix, tmp))
  }
  for(x in 1:ncol(matrix)) {
    diss <- 1-(1.05^(sum(matrix[,x])*-0.02))
    diss <- diss*100
    tmpDiss <- as.matrix(cbind(tmpDiss, diss))
  }
  dissimilarity <- as.matrix(rbind(dissimilarity, tmpDiss))
  tmpDiss <- c()
  matrix <- c()
}

colnames(dissimilarity) <- colnames(normalized)
rownames(dissimilarity) <- colnames(normalized)
write.csv(dissimilarity, file = "dissimilaridade.txt", sep = "\t")

# Hierarchical Clustering and Phylogenetic Tree Generation
hc <- hclust(as.dist(dissimilarity), method = "average")
hc_tree <- as.phylo(hc)
write.tree(hc_tree, file = "M144_dendrogram.newick")

# Plot and Save Dendrogram
pdf(file = "M144.pdf", width = 11, height = 6)
plot(hc, main="M144")
dev.off()

# Define Cluster Membership and Save Results
# Set the number of clusters (k) for hierarchical clustering below
# k = 2
cluster <- rect.hclust(hc, k = 2, border = c("#ff0000", "#4d88ff"))

# Export Cluster Data
cluster_data <- matrix(nrow = 0, ncol = 2)
count <- 1
for (i in cluster) {
  for (x in 1:length(i)) {
    cluster_data <- rbind(cluster_data, c(count, names(i[x])))
  }
  count <- count + 1
}
write.table(cluster_data, file = "clusters_M144.txt", sep = "\t", quote = FALSE, col.names = c("cluster", "indiv"))

# Color Settings for Cluster Plotting
palette <- c("#ff0000", "#4d88ff")
sample_cluster_color <- matrix(nrow = 0, ncol = 4)
count <- 1

for(i in cluster) {
  hex_color <- palette[count]
  for(x in 1:length(cluster[[count]])) {
    sample_cluster_color <- rbind(sample_cluster_color, c(count, names(cluster[[count]])[x], hex_color))
  }
  count <- count + 1  
}

# Initialize and Plot HRM Curve with Clustering Color
plot144 <- plot_ly(x = filter$Temp)
for(i in 1:nrow(sample_cluster_color)) {
  plot144 <- plot144 %>% add_lines(y = filter[[sample_cluster_color[i, 2]]],
                                   color = I(sample_cluster_color[i, 3]),
                                   name = sample_cluster_color[i, 2],
                                   line = list(width = 3))
}

# Customize Layout
plot144 <- plot144 %>%
  layout(
    title = list(text = 'M144/k2', font = list(size = 20)),
    xaxis = list(title = 'Temperature (ºC)', titlefont = list(size = 20), tickfont = list(size = 20)),
    yaxis = list(title = 'Fluorescence (-dF/dT)', titlefont = list(size = 20), tickfont = list(size = 20)),
    legend = list(font = list(size = 14))
  )

orca(plot144, "df_M144.pdf", width=700, height=300)
