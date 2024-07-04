

# """
#
# Filename: SNF Analyses and Independent Replication Validation.R
# Version: 2024.04.2+764
# Author: Sanaz Nazari
# Date Created: July 3, 2024
# Last Modified: July 3, 2024
# Description: This script performs 3 different SNFs on intake, outcome, and mixed data,
#              in addition to the SNF for independent replication validation.
#              The summary of 4 SNFs is provided at the end.
#
# """

#--------------------------
#Similarity Network Fusion

# Required libraries
library(SNFtool)
library(dplyr)

k <- 20        # number of neighbors, usually (10~30)
alpha <- 0.5   # hyperparameter, usually (0.3~0.8)
t <- 20        # Number of Iterations, usually (10~20)
c <- 3         # number of clusters

# Helper function to extract and normalize data
extract_and_normalize <- function(df, cols) {
  temp_df <- data.frame(df[cols])
  normalized_df <- standardNormalization(temp_df)
  normalized_df
}

# Column assignment
select_columns <- function(data, indices) {
  data[indices]
}

# Set dimnames
set_dimnames <- function(df, rows, names) {
  temp_df <- df
  dimnames(temp_df) <- list(1:rows, names)
  temp_df
}

# Create layers
create_layer <- function(data, cols, ss, names) {
  selected_cols <- select_columns(data, cols)
  layer <- extract_and_normalize(selected_cols)
  layer <- set_dimnames(layer, ss, names)
  layer
}

# Layer selection functions
get_fmri_layer <- function(data, ss) {
  create_layer(data, c("col3", "col4", "col5", "col6"), ss, c("var1", "var2", "var3", "var4"))
}

get_clinical_layer <- function(data, cols, names, ss) {
  create_layer(data, cols, ss, names)
}

# Final list creation
combine_layers <- function(layers) {
  intermediate <- list()
  for (i in seq_along(layers)) {
    intermediate[[i]] <- layers[[i]]
  }
  names(intermediate) <- c("a", "b", "c")
  intermediate
}

# Function to read and preprocess data
preprocess_data <- function(file_path, ss, intake = FALSE) {
  data <- read.csv(file_path, header = TRUE, na.strings = c("NULL", "NA", ""))
  
  # Ensure unique column names
  colnames(data)[3:6] <- c("col3", "col4", "col5", "col6")
  colnames(data)[colnames(data) == "rl.1"] <- "col1_rl1"
  colnames(data)[colnames(data) == "el.1"] <- "col2_el1"
  colnames(data)[colnames(data) == "adap.1"] <- "col9_adap1"
  colnames(data)[colnames(data) == "rl"] <- "col7_rl"
  colnames(data)[colnames(data) == "el"] <- "col8_el"
  colnames(data)[colnames(data) == "adap"] <- "col10_adap"
  
  # Layers
  layers <- list()
  layer_funcs <- list(get_fmri_layer, get_clinical_layer, get_clinical_layer)
  clinical_layer_cols <- list(c("col1_rl1", "col2_el1"), "col9_adap1")
  clinical_layer_names <- list(c("var5", "var6"), "var7")
  
  for (i in seq_along(layer_funcs)) {
    if (i == 1) {
      layers[[i]] <- layer_funcs[[i]](data, ss)
    } else {
      cols <- clinical_layer_cols[[i - 1]]
      names <- clinical_layer_names[[i - 1]]
      if (!intake && i > 1) {
        if (i == 2) {
          cols <- c("col7_rl", "col8_el")
          names <- c("var5", "var6")
        } else {
          cols <- "col10_adap"
          names <- "var7"
        }
      }

      layers[[i]] <- layer_funcs[[i]](data, cols, names, ss)
    }
  }
  
  combined_layers <- combine_layers(layers)
  data_layers <- list(data = data)
  data_layers <- c(data_layers, combined_layers)
  data_layers
}

# Function to construct similarity graphs and fuse them
construct_similarity_graphs <- function(layers, k, alpha, t) {
  
  # Helper function to calculate distances
  calculate_distances <- function(df) {
    (dist2(as.matrix(df), as.matrix(df)))^(1/2)
  }
  
  # Helper function to create affinity matrix
  create_affinity_matrix <- function(dist, k, alpha) {
    affinityMatrix(dist, k, alpha)
  }
  
  # Loop through layers
  distance_matrices <- lapply(layers, calculate_distances)
  affinity_matrices <- lapply(distance_matrices, create_affinity_matrix, k = k, alpha = alpha)
  
  # Fuse graphs
  w <- SNF(affinity_matrices, k, t)
  w
}

# Function to summarize data
summarize_data <- function(data, intake = FALSE) {
  
  # Helper function to perform summarization
  perform_summary <- function(df, cols) {
    df %>%
      group_by(labels) %>%
      dplyr::summarise(across(all_of(cols), ~ round(mean(.x, na.rm = TRUE), 2)), n = n())
  }
  
  if (intake) {
    cols <- c("col3", "col4", "col5", "col6", "col1_rl1", "col2_el1", "col9_adap1")
  } else {
    cols <- c("col3", "col4", "col5", "col6", "col7_rl", "col8_el", "col10_adap")
  }
  
  summary <- perform_summary(data, cols)
  summary
}

# Function to process each dataset
process_dataset <- function(file_path, ss, intake = FALSE) {
  preprocessed_data <- preprocess_data(file_path, ss, intake)
  layers <- list(preprocessed_data$a, preprocessed_data$b, preprocessed_data$c)
  w <- construct_similarity_graphs(layers, k, alpha, t)
  labels <- spectralClustering(w, c)
  estimateNumberOfClustersGivenGraph(w, NUMC = 2:10)
  preprocessed_data$data$labels <- labels
  summary <- summarize_data(preprocessed_data$data, intake)
  colnames(summary) <- c("Clusters", "mean.LHfrontal", "mean.RHfrontal", "mean.LHtemporal", "mean.RHtemporal", 
                         "mean.Receptive Language", "mean.Expressive Language", "mean.Adaptive Behavior", "n")
  summary
}

# Processing each dataset
summary_intake <- process_dataset("data-fmri-clin-et-last-1st-81.csv", 81, intake = TRUE)
summary_outcome <- process_dataset("data-fmri-clin-et-last-1st-81.csv", 81)
summary_mixed <- process_dataset("data-fmri-clin-et-last-1st-mixed-137.csv", 137)
summary_replication <- process_dataset("data-nhb-fmri-clin-mixed-last-42.csv", 42)

# Display summaries
# Please note that cluster numbers were created by SNF and the order was changed in the manuscript
print(summary_intake)
print(summary_outcome)
print(summary_mixed)
print(summary_replication)








