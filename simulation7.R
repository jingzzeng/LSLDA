## ---------------------------------------- ##
# This file implements the LSLDA algorithm on simulated data generated from Model (M7) (Section 5).
## ---------------------------------------- ##
rm(list = ls())
# setwd("path/to/LSLDA_demo")   # Set the path to the folder "LSLDA_demo"
source("utility.R")
source("lslda.R")
# #############################################
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

## Model setting
p <- 600
K <- 4
r <- 1
s <- 10
ind <- 1:s
ratio <- rep(1, K)
prior <- ratio/sum(ratio)
N_train <- 30 * ratio
N_test <- 150 * ratio

Sigma <- diag(1, p, p)
Sigma[1:500, 1:500] <- kronecker(diag(10), CS(0.3, 50))
theta <- matrix(0, p, K-1)
theta[1:5,1:3] <- 1
theta[6:10,1:3] <- -1
U <- Sigma %*% theta
Mu <- cbind(rep(0,p), U)
beta <- svd(theta)$u[,1:r, drop = FALSE]
fused_class <- c(2,3,4)

dist <- subspace(svd(Mu)$u[,1:r, drop=FALSE], beta)

## main program
foo <- function(i){
  cat("Time", i, '\n')
  data <- data_gen(N_train, N_test, Mu, Sigma)
  
  x_train <- data$x_train
  y_train <- data$y_train
  x_val <- data$x_val
  y_val <- data$y_val
  x_test <- data$x_test
  y_test <- data$y_test
  
  # ---------------  Bayes error ---------------- #
  fit_bayes <- Bayes_pred(x_test, prior, Mu, Sigma)
  class_bayes <- fit_bayes$class
  err_bayes <-  1 - score_func(class_bayes, y_test, metric = "fused", fused_class = fused_class)
  bayes_result <- err_bayes
  names(bayes_result) <- "err_bayes"
  posterior_bayes <- fit_bayes$posterior
  
  # ---------------- lslda --------------- #
  fit_lslda <- lslda_val(x_train, y_train, x_val, y_val, metric = "loglike", lam1_maxfac=1, n1 = 20, lam2_maxfac = 1)
  beta_lslda <- fit_lslda$beta
  r_lslda <- fit_lslda$rank
  if(is.null(beta_lslda)){
    lslda_result <- rep(NA_real_, 6)
  }
  else{
    ind_hat <- which(apply(beta_lslda, 1, function(x) any(x!=0))) # indices of selected variables
    TPR_lslda <- sum(ind_hat %in% ind)/length(ind)
    FPR_lslda <- ifelse(p == length(ind), 0, sum(!(ind_hat %in% ind))/(p-length(ind)))
    dist_lslda <- subspace(beta, beta_lslda) # subspace distance
    pred_lslda <- predict_lslda(x_train, y_train, beta_lslda, r_lslda, x_test)
    err_lslda <- 1 - score_func(pred_lslda$class, y_test, metric = "fused", fused_class = fused_class)
    kld_lslda <- score_func(pred_lslda$posterior, posterior_bayes, metric = "kld")
    lslda_result <- c(err_lslda, kld_lslda, dist_lslda, TPR_lslda, FPR_lslda, r_lslda)
  }
  names(lslda_result) <- c("err_lslda", "kld_lslda", "dist_lslda", "TPR_lslda", "FPR_lslda",  "r_lslda")
  
  lslda_result
}

times <- 2
output <- lapply(seq_len(times), foo)
output <- do.call(rbind, output)
