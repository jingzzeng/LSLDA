## ---------------------------------------- ##
# This file implements the LSLDA algorithm on simulated data generated from Model (M1) (Section 5).
## ---------------------------------------- ##
rm(list = ls())           # Clear the current environment
# setwd("path/to/LSLDA_demo")   # Set the path to the folder "LSLDA_demo"
source("utility.R")       # Load necessary files for the implementation of SEAS
source("lslda.R")         
# #############################################
RNGkind("L'Ecuyer-CMRG")  # Set the kind of Random Number Generator.
set.seed(1)               # Random seed

## ------------ Model settings ------------ ##
p <- 3000                                    # Dimension of X
K <- 4                                       # The number of class
r <- 2                                       # Discriminant rank
s <- 10                                      # Sparsity level
ind <- 1:s                                   # Active set
ratio <- rep(1, K)                           
prior <- ratio/sum(ratio)                    # Prior
N_train <- 30 * ratio                        # The sample size in each class (training set)
N_test <- 5 * N_train                        # The sample size in each class (test set)
Sigma <- diag(1, p, p)                       # The common covariance
Sigma[1:500, 1:500] <- AR(0.5, 500)
theta <- matrix(0, p, K-1)                   # Theta
theta[1:5,1] <- 0.8
theta[6:10,2] <- 0.8
theta[,3] <- theta[,1] + theta[,2]
U <- Sigma %*% theta
Mu <- cbind(rep(0,p), U)                     # The mean vectors
beta <- svd(theta)$u[,1:r, drop = FALSE]     # The discriminant basis
## ---------------------------------------- ##

## ------------ Main function ------------ ##
foo <- function(i){
  cat("Time", i, '\n')
  data <- data_gen(N_train, N_test, Mu, Sigma)    # Generate data set
  
  x_train <- data$x_train
  y_train <- data$y_train
  x_val <- data$x_val
  y_val <- data$y_val
  x_test <- data$x_test
  y_test <- data$y_test
  
  # ---------------- lslda --------------- #
  fit_lslda <- lslda_val(x_train, y_train, x_val, y_val, metric = "acc")
  beta_lslda <- fit_lslda$beta
  r_lslda <- fit_lslda$rank
  if(is.null(beta_lslda)){
    lslda_result <- rep(NA_real_, 5)
  }
  else{
    ind_hat <- which(apply(beta_lslda, 1, function(x) any(x!=0))) # Estimated active set
    TPR_lslda <- sum(ind_hat %in% ind)/length(ind)
    FPR_lslda <- ifelse(p == length(ind), 0, sum(!(ind_hat %in% ind))/(p-length(ind)))
    dist_lslda <- subspace(beta, beta_lslda) # subspace distance
    pred_lslda <- predict_lslda(x_train, y_train, beta_lslda, r_lslda, x_test) # prediction
    class_lslda <- pred_lslda$class
    err_lslda <- 1- score_func(class_lslda, y_test, metric = "acc")
    lslda_result <- c(err_lslda, dist_lslda, TPR_lslda, FPR_lslda, r_lslda)
  }
  names(lslda_result) <- c("err_lslda", "dist_lslda", "TPR_lslda", "FPR_lslda",  "r_lslda")
  
  lslda_result
}
## ---------------------------------------- ##

times <- 2                                 # Number of replicates
output <- lapply(seq_len(times), foo)      # Implementation of main function
output <- do.call(rbind, output)           # times x 5 output matrix
