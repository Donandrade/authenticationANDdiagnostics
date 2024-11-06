# This script generates a binary matrix from an input file containing data on clusters and individuals.
# The resulting binary matrix can be useful for visualization and annotation in tools such as iTOL.

library(tidyverse)

# Function to generate the binary matrix
generate_binary_matrix <- function(input_file, output_file = NULL) {
  # Read the CSV file
  data <- read.table(input_file, header = TRUE, sep = "\t")
  
  # Check if the column names are correct
  if (!all(c("cluster", "indiv") %in% colnames(data))) {
    stop("The file must contain the columns 'cluster' and 'indiv'.")
  }
  
  # Identify all unique individuals
  unique_individuals <- unique(data$indiv)
  
  # Create the binary matrix
  binary_matrix <- matrix(0, nrow = length(unique(data$cluster)), ncol = length(unique_individuals))
  rownames(binary_matrix) <- unique(data$cluster)
  colnames(binary_matrix) <- unique_individuals
  
  # Populate the binary matrix: 1 indicates the presence of the individual in the corresponding cluster
  for (row in 1:nrow(data)) {
    cluster <- as.character(data$cluster[row])
    individual <- data$indiv[row]
    binary_matrix[cluster, individual] <- 1
  }
  
  # Convert to a data.frame for easier visualization
  binary_df <- as.data.frame(binary_matrix)
  print(binary_df)
  
  # Save the binary matrix to a CSV file if specified
  if (!is.null(output_file)) {
    write.csv(t(binary_df), output_file, row.names = TRUE)
  }
  
  return(binary_df)
}

# Function to process multiple files
process_multiple_files <- function(file_list) {
  # Read the list of files from the input file
  files <- readLines(file_list)
  
  # Process each line in the list (input file)
  for (file in files) {
    input_file <- paste0(file, ".txt")  # Add the ".txt" suffix to the input file name
    output_file <- paste0("binary_matrix_", file, ".csv")  # Define the output file with the prefix "binary_matrix_"
    
    # Call the function to generate and save the binary matrix
    generate_binary_matrix(input_file, output_file)
  }
}

# File with the list of filenames to be processed (one per line)
file_list <- "file_list.txt"  # Replace with the path to your list file

# Call the function to process the files
process_multiple_files(file_list)

