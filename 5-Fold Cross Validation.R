

# """
#
# Filename: 5-Fold Cross Validation.R
# Version: 2024.04.2+764
# Authors: Sanaz Nazari & Javad Zahiri
# Date Created: July 3, 2024
# Last Modified: July 3, 2024
# Description: This script examines the accuracy of SNF by performing 
#              a repeated 5-fold cross validation.
#
# """

#--------------------
# Required libraries
library(caret)
library(SNFtool)

# Function to run k-fold validation analysis
SNF.repeated.cross.validation <- function(no.of.repeat, fold.No, layers.list, cluster.index.vctr) {
  accuracy.vctr.repeated.k.fold.CV <- c()
  
  for (i in 1:no.of.repeat) {
    accuracy.vctr.k.fold.CV <- 
      SNF.cross.validation(fold.No = 5, 
                           layers.list = layers.list, 
                           cluster.index.vctr = cluster.index.vctr)
    
    accuracy.vctr.repeated.k.fold.CV <- 
      c(accuracy.vctr.repeated.k.fold.CV, accuracy.vctr.k.fold.CV)
  }
  return(accuracy.vctr.repeated.k.fold.CV)
}

# Function of cross validation
SNF.cross.validation <- function(fold.No, layers.list, cluster.index.vctr) {
  flds <- createFolds(1:length(cluster.index.vctr), k = fold.No, list = TRUE, returnTrain = FALSE)
  accuracy.vctr <- c()
  
  i <- 1
  for (fold.index in flds) {

    # Prediction of the subtype of new subjects
    testSampleIndexVctr <- fold.index
    
    no.of.train.samples <- length(cluster.index.vctr) - length(testSampleIndexVctr) 
    no.of.test.samples <- length(cluster.index.vctr) - no.of.train.samples
    train <- lapply(layers.list, function(x) x[-testSampleIndexVctr, ]) 
    test <- lapply(layers.list, function(x) x[testSampleIndexVctr, ]) 
    
    train[[3]] <- data.frame(adap = train[[3]])
    test[[3]] <- data.frame(adap = test[[3]])
    
    train <- lapply(seq_along(train), function(i) {
      rownames(train[[i]]) <- 1:no.of.train.samples
      train[[i]]
    })
    
    test <- lapply(seq_along(test), function(i) {
      rownames(test[[i]]) <- 1:no.of.test.samples
      test[[i]]
    })
    
    # labels of the training data
    groups <- cluster.index.vctr[-testSampleIndexVctr] 
    
    # Set the other parameters
    K <- 20
    alpha <- 0.5
    t <- 20
    method <- TRUE
    
    newcluster.index.vctr.from.SNF <- groupPredict(train, test, groups, K, alpha, t, method = 1)
    
    # Compute the prediction accuracy
    accuracy.for.this.fold <- sum(cluster.index.vctr[testSampleIndexVctr] == newcluster.index.vctr.from.SNF[-(1:no.of.train.samples)]) / no.of.test.samples
    accuracy.vctr <- c(accuracy.vctr, accuracy.for.this.fold)
    i <- i + 1
  }
  return(accuracy.vctr)
}


# Function to load data and perform repeated cross-validation
perform_repeated_cross_validation <- function(file_path, ROI_var, clinic_var, label_col, no_of_repeat = 10, fold_no = 5) {
  # Load the data
  SNF_cluster <- read.csv(file_path, header = TRUE, na.strings = c("NULL", "NA", ""))
  
  # Create layers list
  layers.list.outcome <- list(SNF_cluster[, ROI_var], 
                              SNF_cluster[, clinic_var], 
                              data.frame(adap = SNF_cluster[, "adap"]))
  
  # Perform repeated cross-validation
  accuracy.vctr.repeated.k.fold.CV <- SNF.repeated.cross.validation(no.of.repeat = no_of_repeat, 
                                                                    fold.No = fold_no, 
                                                                    layers.list = layers.list.outcome, 
                                                                    cluster.index.vctr = SNF_cluster[, label_col])
  # Return the mean accuracy
  return(mean(accuracy.vctr.repeated.k.fold.CV))
}

# Define variables
ROI_var <- c("LHfrontal", "RHfrontal", "LHtemporal", "RHtemporal")
clinic_var <- c("rl", "el")
label_col <- "labels.last"

# Perform repeated cross-validation for outcome data
outcome_data_path <- "data-fmri-clin-et-last-1st-81.csv"
outcome_accuracy <- perform_repeated_cross_validation(outcome_data_path, ROI_var, clinic_var, label_col)
print(paste("Mean accuracy for Outcome data:", round(outcome_accuracy,2)))

# Perform repeated cross-validation for mixed data
mixed_data_path <- "data-fmri-clin-et-last-mixed-137.csv"
mixed_accuracy <- perform_repeated_cross_validation(mixed_data_path, ROI_var, clinic_var, label_col)
print(paste("Mean accuracy for Mixed data:", round(mixed_accuracy,2)))



