# Analysis Script for Cluster Bootstrap in SNP Data
# This script performs hierarchical clustering with bootstrap resampling 
# on a SNP binary matrix for both rows and columns. The resulting trees 
# are saved in Newick format with bootstrap values, making them compatible 
# with iTOL (Interactive Tree of Life) for further annotation and visualization.

# Install necessary packages (run only once)
# install.packages("ape")
# install.packages("pvclust")

# Load necessary libraries
library(ape)
library(pvclust)

# Load data from CSV file
data <- read.table("fig_S2/data_matrix_binary_Fig_S2.tsv", row.names = 1, header = TRUE, sep = "\t", check.names = TRUE)

# Convert data to matrix
data_matrix <- as.matrix(data)

# Uncomment to remove the first column from the matrix
# data_matrix <- data_matrix[, -1, drop = FALSE]

# Perform bootstrap analysis with pvclust for rows
fit_lines <- pvclust(t(data_matrix), method.hclust = "average", method.dist = "euclidean", nboot = 1000)

# Function to convert hclust to phylo object, including node names. This part of the code came from the StackOverflow topic 
# "How to append bootstrapped values of cluster's (tree) nodes in NEWICK format in R" 
# Link: https://stackoverflow.com/questions/22749634/how-to-append-bootstrapped-values-of-clusters-tree-nodes-in-newick-format-in
as.phylo.hclust.with.nodenames <- function (x, nodenames, ...) {
  N <- dim(x$merge)[1]
  edge <- matrix(0L, 2 * N, 2)
  edge.length <- numeric(2 * N)
  node <- integer(N)
  node[N] <- N + 2L
  cur.nod <- N + 3L
  j <- 1L
  for (i in N:1) {
    edge[j:(j + 1), 1] <- node[i]
    for (l in 1:2) {
      k <- j + l - 1L
      y <- x$merge[i, l]
      if (y > 0) {
        edge[k, 2] <- node[y] <- cur.nod
        cur.nod <- cur.nod + 1L
        edge.length[k] <- x$height[i] - x$height[y]
      } else {
        edge[k, 2] <- -y
        edge.length[k] <- x$height[i]
      }
    }
    j <- j + 2L
  }
  if (is.null(x$labels))
    x$labels <- as.character(1:(N + 1))
  node.lab <- nodenames[order(node)]
  obj <- list(edge = edge, edge.length = edge.length/2, tip.label = x$labels, 
              Nnode = N, node.label = node.lab)
  class(obj) <- "phylo"
  reorder(obj)
}

# Calculate bootstrap values for rows
bootstraps_lines <- (round(fit_lines$edges, 2) * 100)[, 1:2]

# Convert hclust object to phylo with node names for rows
phylo_lines <- as.phylo.hclust.with.nodenames(fit_lines$hclust, nodenames = bootstraps_lines[, 2])

# Save Newick file containing bootstrap values calculated by pvclust for rows
write.tree(phylo_lines, file = "fig_S2/newick_lines_snps_matrix.nwk", tree.names = TRUE, digits = 2)

# Perform bootstrap analysis with pvclust for columns
fit_columns <- pvclust(data_matrix, method.hclust = "average", method.dist = "euclidean", nboot = 1000)

# Display the hclust result for columns
fit_columns$hclust

# Save the result in Newick format for columns
# Extract the resulting tree directly from the fit_columns object and save as Newick
tree <- as.phylo(fit_columns$hclust)

# Display the tree
tree

# Write Newick file for columns
write.tree(tree, file = "fig_S2/newick_columns_snps.nwk", tree.names = TRUE, digits = 2)

