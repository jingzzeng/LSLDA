library(ggplot2)

# Model fitting with lslda function using a separate validation dataset.
lslda_val <- function(x_train, y_train, x_val, y_val, metric = c("accuracy", "gmean", "loglikelihood", "fused"), lam1 = NULL, lam2 = NULL, plot = FALSE, ...){
  metric <- match.arg(metric)
  fit <- lslda(x_train, y_train, lam1, lam2, ...)
  pred <- predict.lslda(fit, x_val)
  rank_l <- fit$rank
  s_l <- fit$s
  beta_l <- fit$beta
  step_l <- fit$step
  time_l <- fit$time
  lam1 <- fit$lam1
  lam2 <- fit$lam2
  ##
  ind_sparse <- which((s_l >= 5) & (s_l <= 15))
  time_sparse_l <- time_l[ind_sparse]
  step_sparse_l <- step_l[ind_sparse]
  s_sparse_l <- s_l[ind_sparse]
  ##
  
  if(metric == "loglikelihood"){
    score <- sapply(pred$posterior, score_func, y_val, metric = metric)
  }else{
    score <- apply(pred$class, 2, score_func, y_val, metric = metric)
  }
  
  id <- which(score == max(score))
  if(all(is.na(rank_l[id]))){
    id <- id[1]
  }else if(any(is.na(rank_l[id]))){
    id <- id[!is.na(rank_l[id])]
    id <- id[rank_l[id] == min(rank_l[id], na.rm = TRUE)]
    id <- id[which.min(s_l[id])]
  }else{
    id <- id[rank_l[id] == min(rank_l[id], na.rm = TRUE)]
    id <- id[which.min(s_l[id])]
  }
  
  if(plot){
    dat <- data.frame(x = 1:length(score), y = score, rank = as.factor(rank_l))
    g <- ggplot(dat, aes(x = x, y = y, col = rank))+
      geom_point(size = 1)+
      labs(x="", y="Validation score")+
      theme_bw()
    print(g)
  }
  
  n1 <- length(lam1)
  n2 <- length(lam2)
  id_lam1 <- ceiling(id/(n2))
  id_lam2 <- id-(id_lam1-1)*n2
  lam1_max <- lam1[id_lam1]
  lam2_max <- lam2[id_lam2]
  score_max <- score[id]
  rank <- rank_l[id]
  s <- s_l[id]
  beta <- beta_l[[id]]
  
  output <- list(beta = beta, rank = rank, s = s, score = score_max, lam1 = lam1, lam2 = lam2, lam1_max = lam1_max, lam2_max = lam2_max, step = step_l, time = time_l, s_sparse = s_sparse_l, step_sparse = step_sparse_l, time_sparse = time_sparse_l)
}

# Implementing LSLDA with tuning parameters sequencees
lslda <- function(x, y, lam1 = NULL, lam2 = NULL, ...){
  sig_mu <- MU(x, y)
  M <- sig_mu$M
  U <- sig_mu$U
  largesteig <- sig_mu$largesteig
  Mhalf <- sig_mu$Mhalf
  fit <- splitting_algo(M, U, lam1, lam2, largesteig =  largesteig, Mhalf = Mhalf, ...)
  output <- list(x = x, y = y, beta = fit$beta_l, rank = fit$rank_l, s = fit$s_l, lam1 = fit$lam1_l, lam2 = fit$lam2_l, step = fit$step_l, time = fit$time_l)
  output
}

# Cross-validation function for LSLDA
cv.lslda <- function(x, y, nfolds = 5, stratified = FALSE, metric = c("accuracy", "gmean", "loglikelihood", "fused"), lam1 = NULL, lam2 = NULL, ...){
  metric <- match.arg(metric)
  x_old <- x; y_old <- y
  if(is.factor(y)){y <- as.numeric(y)}
  y <- drop(y)
  n <- nrow(x)
  p <- ncol(x)
  ord <- order(y)
  y <- y[ord]
  x <- x[ord,]
  fit <- lslda(x, y, lam1, lam2, ...)
  lam1 <- fit$lam1
  lam2 <- fit$lam2
  n1 <- length(lam1)
  n2 <- length(lam2)
  if(stratified){
    cls_ss <- table(y)
    foldid <- c()
    for(i in cls_ss){
      foldid <- c(foldid, sample(rep(sample(seq_len(nfolds)), length = i)))
    }
  }else{
    foldid <- sample(rep(sample(seq_len(nfolds)), length = n))
  }
  cv_result <- matrix(0, nfolds, n1*n2)
  for (i in seq_len(nfolds)) {
    cat("Fold", i, "\n")
    which <- foldid == i
    x_train <- x[!which, , drop = FALSE]
    y_train <- y[!which]
    x_val <- x[which, , drop = FALSE]
    y_val <- y[which]
    fit_cv <- lslda(x_train, y_train, lam1, lam2, ...)
    pred <- predict.lslda(fit_cv, x_val)
    if(metric == "loglikelihood"){
      score <- sapply(pred$posterior, score_func, y_val, metric = metric)
    }else{
      score <- apply(pred$class, 2, score_func, y_val, metric = metric)
    }
    cv_result[i,] <- score
  }
  cvm <- colMeans(cv_result, na.rm = TRUE)
  rank_l <- fit$rank
  s_l <- fit$s
  beta_l <- fit$beta
  id <- which(cvm == max(cvm))
  if(all(is.na(rank_l[id]))){
    id <- id[1]
  }else if(any(is.na(rank_l[id]))){
    id <- id[!is.na(rank_l[id])]
    id <- id[rank_l[id] == min(rank_l[id], na.rm = TRUE)]
    id <- id[which.min(s_l[id])]
  }else{
    id <- id[rank_l[id] == min(rank_l[id], na.rm = TRUE)]
    id <- id[which.min(s_l[id])]
  }
  rank <- rank_l[id]
  s <- s_l[id]
  beta <- beta_l[[id]]
  output <- list(x = x_old, y = y_old, beta = beta, rank = rank, s = s, cvm = cvm, lam1 = lam1, lam2 = lam2, lslda.fit = fit)
  output
}

# The core function implementing the LSLDA algorithm
splitting_algo <- function(M, U, lam1_l, lam2_l, largesteig = NULL, Mhalf = NULL, ...){
  if(is.null(M)){
    if(is.null(Mhalf)) stop("M and Mhalf are both missing.")
  }
  opts <- list(...)
  if(is.null(lam1_l)){
    if(is.null(opts$n1)) opts$n1 <- 10
    if(is.null(opts$lam1_maxfac)) opts$lam1_maxfac <- 0.5
    if(is.null(opts$lam1_minfac)) opts$lam1_minfac <- 0.1
    n1 <- opts$n1
    lam1_maxfac <- opts$lam1_maxfac
    lam1_minfac <- opts$lam1_minfac
    Urow_norm <- apply(U, 1, function(x){sqrt(sum(x^2))})
    lam1_max <- max(Urow_norm)*lam1_maxfac
    # lam1_min <- max(Urow_norm)*lam1_minfac
    lam1_min <- lam1_max*lam1_minfac
    lam1_l <- seq(lam1_max, lam1_min, length.out = n1)
  }

  if(is.null(lam2_l)){
    if(is.null(opts$n2)) opts$n2 <- 15
    if(is.null(opts$lam2_maxfac)) opts$lam2_maxfac <- 0.5
    if(is.null(opts$lam2_minfac)) opts$lam2_minfac <- 0.1
    n2 <- opts$n2
    lam2_maxfac <- opts$lam2_maxfac
    lam2_minfac <- opts$lam2_minfac
    U_sv <-  svd(U)$d
    lam2_max <- U_sv[1]*lam2_maxfac
    # lam2_min <- U_sv[1]*lam2_minfac
    lam2_min <- lam2_max*lam2_minfac
    lam2_l <- seq(lam2_max, lam2_min, length.out = n2)
  }

  if(is.null(opts$pmax_fac)) opts$pmax_fac <- 0.5
  if(is.null(opts$tol)) opts$tol <- 1e-3
  if(is.null(opts$maxiter)) opts$maxiter <- 1e3
  if(is.null(opts$r_thrd)) opts$r_thrd <- 1e-3

  if(is.null(largesteig)){
    gamma <- 1.99/(eigen(M, symmetric = TRUE)$values[1])
  }else{
    gamma <- 1.99/largesteig
  }

  pmax_fac <- opts$pmax_fac
  tol <- opts$tol # covergence tolerance
  maxiter <- opts$maxiter # maximal iteration
  r_thrd <- opts$r_thrd # rank threshold
  p <- NROW(U)
  m <- NCOL(U)
  lambda <- 1
  n1 <- length(lam1_l)
  n2 <- length(lam2_l)
  if(!is.null(Mhalf)) Mhalft <- t(Mhalf)

  rank_l <- rep(NA_integer_, n1*n2)
  beta_l <- vector("list", n1*n2)
  iter_l <- rep(NA_integer_, n1*n2)
  s_l <- rep(NA_integer_, n1*n2)
  sv_l <- vector("list", n1*n2)
  step_l <- rep(NA_real_, n1*n2)
  time_l <- rep(NA_real_, n1*n2)
  for(i in seq_along(lam1_l)){
    lam1 <- lam1_l[i]
    for(j in seq_along(lam2_l)){
      lam2 <- lam2_l[j]
      z <- matrix(0, p, m)
      t <- 0
      # main iteration
      start_time <- Sys.time() 
      step <- 0
      code <- 0
      if(is.null(Mhalf)){
        while(t < maxiter){
          step <- step + 1
          x_f_out <- prox_gs(z, gamma * lam1) # step 1: group sparsity
          x_f <- x_f_out$B_out
          nz_ind <- x_f_out$nz_ind
          if(is.null(nz_ind)){
            tmp <- 2*x_f - z + gamma * U
          }else{
            tmp <- 2*x_f - z - gamma * (M[, nz_ind, drop=FALSE] %*% x_f[nz_ind, ,drop=FALSE] - U)
          }
          x_g <- prox_lr(tmp, gamma * lam2) # step 2: low-rank
          z <- z + lambda * (x_g - x_f) # step 3: update B
          t <- t + 1
          change <- sqrt(sum((x_g - x_f)^2))/(1 + sqrt(sum(z^2)))
          if(is.infinite(change)){code <- -1; break}
          else if(change < tol){code <- 1; break} # convergence criterion: relative change.
        }
      }
      else{
        while(t < maxiter){
          step <- step + 1
          x_f_out <- prox_gs(z, gamma * lam1) # step 1: group sparsity
          x_f <- x_f_out$B_out
          nz_ind <- x_f_out$nz_ind
          if(is.null(nz_ind)){
            tmp <- 2*x_f - z + gamma * U
          }else{
            tmp <- 2*x_f - z - gamma * (Mhalft %*% (Mhalf[, nz_ind, drop=FALSE] %*% x_f[nz_ind, ,drop=FALSE]) - U)
          }
          x_g <- prox_lr(tmp, gamma * lam2) # step 2: low-rank
          z <- z + lambda * (x_g - x_f) # step 3: update B
          t <- t + 1
          change <- sqrt(sum((x_g - x_f)^2))/(1 + sqrt(sum(z^2)))
          if(is.infinite(change)){code <- -1; break}
          else if(change < tol){code <- 1; break} # convergence criterion: relative change.
        }
      }
      end_time <- Sys.time()
      Bhat <- x_f # final iterate

      if(code == -1) next
      # record
      if((code == 1) || (code == 0)){
        id <- (i-1)*n2 + j
        iter_l[id] <- t
        sv_l[[id]] <- svd(Bhat)$d
        rank <- rank_func(Bhat, r_thrd)
        rank_l[id] <- rank
        step_l[id] <- step
        time_l[id] <- difftime(end_time, start_time, units = "secs")
        if(rank == 0){
          beta <- matrix(0, nrow(Bhat), 1)
        }else{
          beta <- svd(Bhat)$u[,1:rank, drop = FALSE]
        }
        var_ind <- apply(beta, 1, function(x){any(x!=0)})
        s_l[id] <- sum(var_ind)
        beta_l[[id]] <- beta
      }
      if(s_l[id] > p*pmax_fac) {code <- -2; break} # If too dense, break
    }
  }
  output <- list(lam1_l = lam1_l, lam2_l = lam2_l, gamma = gamma, beta_l = beta_l, rank_l = rank_l, s_l = s_l, sv_l = sv_l, step_l = step_l, time_l = time_l)
  output
}

# Prediction function for LSLDA: input the LSLDA object returned from "lslda()" function.
predict.lslda <- function(object, newx){
  x <- object$x
  y <- object$y
  if(is.factor(y)){
    class <- levels(y)
    y <- as.numeric(y)
    out <- predict_lslda(x, y, object$beta, object$rank, newx)
    out$class <- class[out$class]
  }else if(is.numeric(y)){
    out <- predict_lslda(object$x, object$y, object$beta, object$rank, newx)
  }
  list(posterior = out$posterior, class = out$class)
}

# Another prediction function for LSLDA: input x, y, the estimated discriminant basis, the discriminant rank, and the new x. Returns predicted classes and posteriors for the new x.
predict_lslda <- function(x, y, Beta, rank, newx){
  uclass <- sort(unique(y))
  prior <- sapply(uclass, function(t) mean(y==t))
  if(is.atomic(Beta)){
    if(is.null(Beta) || is.na(Beta) || rank == 0){
      # zero matrix or the matrix that is not converged.
      class <- rep(which.max(prior), NROW(newx))
      tmp <- matrix(0, NROW(newx), length(uclass))
      tmp[, which.max(prior)] <- 1
      posterior <- tmp
    }else{
      out <- lda_pred(x, y, Beta, newx)
      class <- out$class
      posterior <- out$posterior
    }
  }
  else if(is.list(Beta)){
    class <- matrix(0, nrow(newx), length(Beta))
    posterior <- vector("list", length(Beta))
    for(i in seq_len(length(Beta))){
      beta <- Beta[[i]]
      r <- rank[[i]]
      if(is.null(beta) || is.na(Beta) || r == 0){
        class[,i] <- which.max(prior)
        tmp <- matrix(0, nrow(newx), length(uclass))
        tmp[, which.max(prior)] <- 1
        posterior[[i]] <- tmp
      }else{
        out <- lda_pred(x, y, beta, newx)
        class[,i] <- out$class
        posterior[[i]] <- out$posterior
      }
    }
  }
  list(class = class, posterior = posterior)
}

# Classical LDA prediction.
lda_pred <- function(x, y, B, newx){
  cls <- unique(y)
  prior <- rep(0, length(cls))
  for (i in seq_along(cls)){
    prior[i] <- mean(y == cls[i])
  }
  proj_x <- x %*% B
  proj_newx <- newx %*% B
  out <- tryCatch({
    fit <- lda(proj_x, y)
    predict(fit, proj_newx)
  },
  error = function(e){
    class <- rep(which.max(prior), NROW(newx))
    posterior <- matrix(0, NROW(newx), length(cls))
    posterior[,which.max(prior)] <- 1
    list(class = class, posterior = posterior)
  })
  list(class = out$class, posterior = out$posterior)
}

# The proximity operator of group sparse penalty.
prox_gs <- function(B, a){
  ones <- matrix(1, nrow = 1, ncol = NCOL(B))
  Brow_norm <- sqrt(tcrossprod(B^2, ones))
  if(a == 0){
    nz_ind <- which(Brow_norm != 0)
    return(list(B_out = B, nz_ind = nz_ind))
  }
  b <- pmax(0, 1-a/Brow_norm)
  if(all(b == 0)){
    nz_ind <- NULL
  }else{
    nz_ind <- which(b != 0)
  }
  B_out <- b * B
  list(B_out = B_out, nz_ind = nz_ind)
}

# The proximity operator of nuclear norm penalty.
prox_lr <- function(B, a){
  svd_out <- svd(B)
  u <- svd_out$u
  v <- svd_out$v
  d <- svd_out$d
  c <- pmax(0, d-a)
  B_out <- u %*% (c * t(v))
  B_out
}

# The estimation function computing the common covariance and the U matrix.
MU <- function(x, y){
  x <- as.matrix(x)
  y <- sort(drop(y))
  class <- unique(y)
  K <- length(class)
  prior <- rep(0, K)
  sigma <- NULL
  for (k in 1:K) {
    prior[k] <- mean(y == class[k])
  }
  p <- as.integer(ncol(x))
  n <- nrow(x)
  ny <- length(y)
  if (ny != n) 
    stop("x and y have different number of observations")
  U <- matrix(0, p, K)
  for (i in 1:K) {
    index <- which(y == class[i])
    n_i <- length(index)
    mu_i <- colMeans(x[index,,drop=FALSE])
    U[, i] <- sqrt(prior[i]) * (mu_i - colMeans(x))
    x[index, ] <- x[index, ] - matrix(mu_i, n_i, p, byrow = TRUE)
  }
  sigma <- NULL
  if(p <= 1000) sigma <- crossprod(x/sqrt(n))
  Mhalf <- x/sqrt(n)
  largesteig <- (svd(x)$d[1])^2/(n)
  outlist <- list(M = sigma, U = U, prior = prior, largesteig = largesteig, Mhalf = Mhalf)
  outlist
}
