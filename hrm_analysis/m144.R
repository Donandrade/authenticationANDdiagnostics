# Install required packages
install.packages("tidyr")
if (!require("processx")) install.packages("processx")
install.packages("MBmca")

# Load necessary libraries
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

###########################################################################
###########################################################################
###########################################################################
###########################################################################
# Read data from file
data <- read.csv("../data/M144_Raw_Data.txt", header = TRUE, sep = "\t", dec = ".",
                 check.names = FALSE)

data2 <- c()

# Filter rows based on the first column's value
for(i in 1:nrow(data)){
  if(data[i,1] > 77 & data[i,1] < 84){
    tmp <- c(data[i,])
    data2 <- as.data.frame(rbind(as.data.frame(data2), as.data.frame(tmp)))
  }
}

# Adjust column names if they start with 'X'
names(data2)[startsWith(names(data2), "X")] <- substring(names(data2)[startsWith(names(data2), "X")], 2)

data <- data2

# Apply meltcurve function on the selected columns
x <- meltcurve(data %>% select(1, 2), norm = TRUE,  span.smooth = 0.05, 
               is.deriv = TRUE, peaklines = TRUE,
               calc.Area = TRUE)

c = x[[1]][1]

c = cbind(c, x[[1]][3])

d = x[[1]][1]
d = cbind(d, x[[1]][["Fluo"]])

# Apply meltcurve on the rest of the columns
for(i in 3:ncol(data)){
  x <- meltcurve(data %>% select(1, i), norm = TRUE,  span.smooth = 0.05, 
                 is.deriv = TRUE, peaklines = TRUE,
                 calc.Area = TRUE)
  c = cbind(c, x[[1]][3])
  d = cbind(d, x[[1]][["Fluo"]])
  
  print (i)
}

colnames(c) = colnames(data)
colnames(d) = colnames(data)

filter <- c()

# Filter rows for specific condition
for(i in 1:nrow(c)){
  if(c[i,1] > 77 & c[i,1] < 84){
    tmp <- c(c[i,])
    print(tmp)
    filter <- as.data.frame(rbind(as.data.frame(filter), as.data.frame(tmp)))
  }
}

colnames(d) <- colnames(data)

normalized = d[,2:ncol(d)]

####################################################
####################################################
####################################################
# Calculate dissimilarity matrix

matrix <- c()
tmpDiss <- c()
dissimilarity <- c()

# Calculate dissimilarity between columns
for(a in 1:ncol(normalized)){
  for(j in 1:ncol(normalized)){
    tmp <- as.matrix((normalized[,a]-normalized[,j])^2)
    matrix <- as.matrix(cbind(matrix, tmp))
  }
  for(x in 1:ncol(matrix)){
    diss <- 1-(1.05^(sum(matrix[,x])*-0.02))
    diss <- diss*100
    tmpDiss <- as.matrix(cbind(tmpDiss, diss))
  }
  dissimilarity <- as.matrix(rbind(dissimilarity, tmpDiss))
  tmpDiss <- c()
  matrix <- c()
}

# Assign column and row names to dissimilarity matrix
colnames(dissimilarity) <- colnames(normalized)
rownames(dissimilarity) <- colnames(normalized)

# Write dissimilarity matrix to a file
write.csv(dissimilarity, file = "dissimilaridade.txt", sep = "\t")

# Perform hierarchical clustering
hc <- hclust(as.dist(dissimilarity), method = "average")

# Save hierarchical clustering plot as PDF
pdf(file = "M144.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6) # The height of the plot in inches

plot(hc, main="M144")

########################################
#### Define the cluster number
########################################

# Rectangular cut for clustering
cluster <- rect.hclust(hc, k = 2, border = 2:5)

dev.off()

# Create an empty matrix to store cluster data
cluster_data <- matrix(nrow = 0, ncol = 2)

count =1

# Populate cluster data matrix with cluster assignments
for (i in cluster) {
  for (x in 1:length(i)) {
    cluster_data <- rbind(cluster_data, c(count, names(i[x])))  # Add row to matrix
  }
  count <- count + 1
}

# Display structure of cluster
str(cluster)

# Write cluster data to a file
write.table(cluster_data, file = "clusters_M144.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,
            col.names =  c("cluster", "indiv"))

#################################
##### Color plant build option. Coloring based on cluster
#################################

# Load RColorBrewer library
library(RColorBrewer)

# Define the number of colors based on number of clusters
num_clusters <- length(cluster)

# Define color palette
palette <- c( "#ff0000", "#4d88ff")

# Initialize empty matrix for sample cluster color assignments
sample_cluster_color <- matrix(nrow = 0, ncol = 4)

# Assign colors to clusters
count <- 1
for(i in cluster){
  hex_color <- palette[count]
  
  for(x in 1:length(cluster[[count]])){
    sample_cluster_color <- rbind(sample_cluster_color,c(count, names(cluster[[count]])[x], hex_color))
  }
  count <- count + 1  
}

# Display sample cluster color assignments
sample_cluster_color

####################################################
####################################################
####################################################
# Set row and column names for filtered data
colnames(data)
rownames(filter) <- filter[,1]
filter
colnames(filter) = colnames(data)

# Initialize the plot
plot <- plot_ly(x = filter$Temp)

# Add lines to the plot for each cluster
for(i in 1:nrow(sample_cluster_color)){
  plot <- plot %>% add_lines(y = filter[[sample_cluster_color[i, 2]]],
                             color = I(sample_cluster_color[i, 3]),
                             name = sample_cluster_color[i, 2],
                             line = list(width = 3))
}

plot

# Adjust layout to increase font sizes
plot <- plot %>%
  layout(
    title = list(text = 'M144/k2', font = list(size = 20)),
    xaxis = list(title = 'Temperature (ºC)', titlefont = list(size = 20), tickfont = list(size = 20)),
    yaxis = list(title = 'Fluorescence (-dF/dT)', titlefont = list(size = 20), tickfont = list(size = 20)),
    legend = list(font = list(size = 14))
  )

# Display the plot
plot

# Save the plot as a PDF
orca(plot, "df_M144.pdf", width=700, height=300)

