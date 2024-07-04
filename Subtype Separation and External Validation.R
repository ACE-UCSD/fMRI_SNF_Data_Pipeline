

# """
#
# Filename: Subtype Separation and External Validation.R
# Version: 2024.04.2+764
# Author: Sanaz Nazari
# Date Created: July 3, 2024
# Last Modified: July 3, 2024
# Description: This script examines subtype separation and external validation by conducting one-way ANOVA 
#              and pairwise comparisons.
#
# """

#--------------------
# Required libraries
library(ez)
library(dplyr)

# Function to check normality using Shapiro-Wilk test
check_normality <- function(data, var, labels) {
  pvalues <- sapply(labels, function(label) {
    shapiro.test(data[[var]][data$labels.last == label])$p.value
  })
  all(pvalues > 0.05)
}

# Function for pairwise comparisons
pairwise_comparisons <- function(data, var, method, labels) {
  if (length(labels) < 2) return(NULL)
  combs <- combn(labels, 2, simplify = FALSE)
  pvalues <- sapply(combs, function(lbl) {
    group1 <- data[[var]][data$labels.last == lbl[1]]
    group2 <- data[[var]][data$labels.last == lbl[2]]
    if (method == "t.test") {
      t.test(group1, group2, var.equal = FALSE)$p.value
    } else {
      kruskal.test(list(group1, group2))$p.value
    }
  })
  p.adjust(pvalues, method = "fdr")
}

# Function to perform ezANOVA and collect results
perform_anova <- function(data, vars) {
  res.anova <- NULL
  for (var in vars) {
    data_copy <- data
    colnames(data_copy)[which(names(data_copy) == var)] <- "var"
    anv <- ez::ezANOVA(data=data_copy, dv=.(var), wid=.(id), type=3, between=.(labels.last), detailed=TRUE, return_aov=TRUE)
    temp <- cbind(var, round(anv$ANOVA$ges, 3), anv$ANOVA$`p<.05`, anv$ANOVA$Effect, round(anv$ANOVA$F, 2),
                  anv$ANOVA$DFn, anv$ANOVA$DFd, round(anv$ANOVA$p, 3))
    res.anova <- rbind(res.anova, temp)
    colnames(data_copy)[which(names(data_copy) == "var")] <- var # Restore original column name
  }
  return(res.anova)
}

# Function to perform full analysis (ANOVA + Pairwise comparisons)
perform_full_analysis <- function(file_path, variables) {
  # Load data
  data <- read.csv(file_path, header=TRUE, na.strings=c("NULL", "NA",""))
  data.anova <- data.frame(data[,c(3:6, 17, 20:21, 15, 22, 12, 43)])
  
  # Convert necessary columns to factors
  data.anova$id <- as.factor(seq(1, nrow(data.anova)))
  data.anova$labels.last <- as.factor(data.anova$labels.last)
  
  # Checking assumptions for ANOVA:
  labels <- unique(as.character(data.anova$labels.last))
  
  # Perform ANOVA on SNF and two external variables
  snf_vars <- variables[1:9]
  res_anova_snf <- perform_anova(data.anova, snf_vars)
  colnames(res_anova_snf) <- c("Variable", "GES",	"p<.05", "Effect", "F-ratio",	"DFn", "DFd" ,	"P-value")
  print(res_anova_snf)
  
  # Perform ANOVA on Eye Tracking variable
  et_vars <- variables[10]
  data.et <- subset(data.anova, !is.na(perc.soc))
  res_anova_et <- perform_anova(data.et, et_vars)
  colnames(res_anova_et) <- c("Variable", "GES",	"p<.05", "Effect", "F-ratio",	"DFn", "DFd" ,	"P-value")
  print(res_anova_et)
  
  # Pairwise comparisons
  results <- list()
  
  for (var in variables) {
    method <- if (check_normality(data.anova, var, labels)) "t.test" else "kruskal"
    adj.p <- pairwise_comparisons(data.anova, var, method, labels)
    if (!is.null(adj.p)) {
      comparisons <- combn(labels, 2, paste, collapse = " vs ")
      results[[var]] <- data.frame(Comparison = comparisons, Adjusted.P = round(adj.p, 3))
    }
  }
  
  # Combine results into a single data frame
  final_results <- do.call(rbind, lapply(names(results), function(var) {
    cbind(Variable = var, results[[var]])
  }))
  
  print(final_results)
  
  return(list(anova_snf = res_anova_snf, anova_et = res_anova_et, pairwise = final_results))
}

# Variables to analyze
variables <- c("LHfrontal", "RHfrontal", "LHtemporal", "RHtemporal", 
               "adap", "rl", "el", "tot", "elc", "perc.soc")

# Perform analysis on outcome data
print("Outcome Data Analysis:")
outcome_results <- perform_full_analysis("data-fmri-clin-et-last-1st-81.csv", variables)

# Perform analysis on mixed data
print("Mixed Data Analysis:")
mixed_results <- perform_full_analysis("data-fmri-clin-et-last-1st-mixed-137.csv", variables)





