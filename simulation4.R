## ---------------------------------------- ##
# This file implements the LSLDA algorithm on simulated data generated from Model (M4) (Section 5).
## ---------------------------------------- ##
rm(list = ls())
# setwd("path/to/LSLDA_demo")   # Set the path to the folder "LSLDA_demo"
source("utility.R")
source("lslda.R")
# #############################################
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# ## Model setting
p <- 3000
K <- 7
r <- 2
s <- 10
ind <- 1:s
ratio <- rep(1, K)
prior <- ratio/sum(ratio)
N_train <- 30 * ratio
N_test <- 5 * N_train

Sigma <- diag(1, p, p)
Sigma[1:500, 1:500] <- AR(0.5, 500)
theta <- matrix(0, p, K-1)
theta[seq(1,s,2),1] <- 2
theta[seq(2,s,2),2] <- -4
tmp <- matrix((1:(K-3))/2, 2, K-3, byrow = TRUE)
beta2 <- theta[,1:2] %*% tmp
theta[,3:(K-1)] <- beta2
U <- Sigma %*% theta
Mu <- cbind(rep(0,p), U)
beta <- svd(theta)$u[,1:r, drop = FALSE]

dist <- subspace(svd(Mu)$u[,1:r], beta)

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
  
  # ---------------- lslda --------------- #
  fit_lslda <- lslda_val(x_train, y_train, x_val, y_val, lam1 = seq(0.5,1,length.out = 10), lam2 = seq(0.3,1.5,length.out = 20), metric = "acc")
  beta_lslda <- fit_lslda$beta
  r_lslda <- fit_lslda$rank
  if(is.null(beta_lslda)){
    lslda_result <- rep(NA_real_, 5)
  }
  else{
    ind_hat <- which(apply(beta_lslda, 1, function(x) any(x!=0))) # indices of selected variables
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

times <- 2
output <- lapply(seq_len(times), foo)
output <- do.call(rbind, output)
