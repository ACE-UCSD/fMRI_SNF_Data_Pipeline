

# """
#
# Filename: Normalized Mutual Information.R
# Version: 2024.04.2+764
# Authors: Sanaz Nazari
# Date Created: July 3, 2024
# Last Modified: July 3, 2024
# Description: This script calculates NMI to examine the robustness of SNF to random removals.
#
# """

#-------------------
# Required libraries
library(SNFtool)

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

# Function to get labels after SNF
get_labels_after_snf <- function(file_path, ss, intake = FALSE) {
  preprocessed_data <- preprocess_data(file_path, ss, intake)
  layers <- list(preprocessed_data$a, preprocessed_data$b, preprocessed_data$c)
  w <- construct_similarity_graphs(layers, k, alpha, t)
  labels <- spectralClustering(w, c)
  labels
}

# Function to perform 100 replications of random % removals and calculate NMI
calculate_nmi <- function(file_path, ss, c, random_ratios, output_file) {
  data <- read.csv(file_path, header = TRUE, na.strings = c("NULL", "NA", ""))
  
  nmi <- list()
  
  for (r in random_ratios) {
    nmi.rr <- list()
    for (i in 1:100) {
      n <- floor(r * nrow(data))
      ind <- sample(seq_len(nrow(data)), size = n)
      data.rr <- data[ind, ]
      rownames(data.rr) <- NULL
      write.csv(data.rr, "data.rr.csv", row.names = FALSE)
      
      #labels <- get_labels_after_snf(file_path, ss, intake = FALSE)
      labels.rr <- get_labels_after_snf("data.rr.csv", n, intake = FALSE)
      data.rr <- data.frame(data.rr, labels.rr)
      
      temp <- calNMI(data.rr$labels.last, data.rr$labels.rr)
      nmi.rr <- append(nmi.rr, list(temp))
      
      print(paste("i =", i, ", r =", r))
    }
    
    temp2 <- nmi.rr
    nmi <- qpcR:::cbind.na(nmi, temp2)
  }
  nmi <- nmi[,-1]
  colnames(nmi) <- c("5%","10%","20%","30%","40%","50%")
  print(nmi)
  
  write.csv(nmi, file = output_file, row.names = FALSE)
}  


# Calculate NMI for 100 replications of random % removals for outcome data
calculate_nmi("data-fmri-clin-et-last-1st-81.csv", 81, c, c(0.95, 0.90, 0.80, 0.70, 0.60, 0.50), "nmi-last-81.csv")

# Calculate NMI for 100 replications of random % removals for mixed data
calculate_nmi("data-fmri-clin-et-last-mixed-137.csv", 137, c, c(0.95, 0.90, 0.80, 0.70, 0.60, 0.50), "nmi-mixed-137.csv")


