library(MASS)

# AR function
AR <- function(rho, p){
  m <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      m[i,j] <- rho**(abs(i-j))
    }
  }
  return(m)
}

# CS function
CS <- function(rho, p){
  m <- matrix(rho,p,p)
  diag(m) <- 1
  return(m)
}

# rank function
rank_func <- function(B, thrd = 1e-3){
  if(is.null(B)) return(NA)
  d <- svd(B)$d
  r <- sum(d >= thrd)
  r
}

# subspace function for matrices with the same column dimension
subspace <- function(A,B){
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*d))
}

# data generator
data_gen <- function(N_train, N_test, Mu, Sigma){
  y_train <- lapply(seq_along(N_train), function(i){rep(i, N_train[i])})
  y_train <- do.call(c, y_train)
  x_train <- x_gen(y_train, Mu, Sigma)
  y_val <- lapply(seq_along(N_train), function(i){rep(i, N_train[i])})
  y_val <- do.call(c, y_val)
  x_val <- x_gen(y_val, Mu, Sigma)
  y_test <- lapply(seq_along(N_test), function(i){rep(i, N_test[i])})
  y_test <- do.call(c, y_test)
  x_test <- x_gen(y_test, Mu, Sigma)
  list(x_train = x_train, y_train = y_train, x_val = x_val, y_val = y_val, x_test = x_test, y_test = y_test)
}

# Generate predictors from conditional normal distribution.
x_gen <- function(y, mu, sigma){
  y <- sort(y)
  n <- length(y)
  p <- NROW(sigma)
  class <- unique(y)
  tab <- table(y)
  mat <- lapply(seq_along(class), function(i){mvrnorm(tab[i], mu[,i], sigma)})
  mat <- do.call(rbind, mat)
  mat
}

# Score function
score_func <- function(pred, y, metric = c("accuracy", "gmean", "loglikelihood", "fused_accuracy", "kld"), fused_class = NULL){
  metric <- match.arg(metric)
  if(metric == "accuracy"){
    score <- mean(pred == y)
  }else if(metric == "gmean"){
    class <- unique(y)
    acc_l <- sapply(class, function(i){sum((y==i) & (pred==i))/sum(y==i)})
    score <- prod(acc_l)^{1/length(class)}
  }else if(metric == "loglikelihood"){
    loglh <- sapply(seq_len(length(y)), function(i){log(pred[i, y[i]])})
    score <- sum(loglh)
  }else if(metric == "fused_accuracy"){
    if(is.list(fused_class)){
      for(classes in fused_class){
        pred[pred %in% classes] <- min(classes)
        y[y %in% classes] <- min(classes)
      }
    }else{
      pred[pred %in% fused_class] <- min(fused_class)
      y[y %in% fused_class] <- min(fused_class)
    }
    score <- mean(pred == y)
  }else if(metric == "kld"){
    if(any(pred == 0)) pred[pred == 0] <- pred[pred == 0] + .Machine$double.eps # avoid Inf
    score <- sum(log(pred/y) * pred)/NROW(pred)
  }
  score
}

# Bayes classifier
Bayes_pred <- function(data, prior, mu, sigma){
  sigma_halfinv <- pracma::sqrtm(sigma)$Binv
  posterior <- sapply(1:NROW(data), function(i){
    rw <- data[i,]
    tmp <-  sigma_halfinv %*% (as.matrix(rw)[,rep(1,K)] - mu)
    tmp2 <- -0.5 * apply(tmp, 2, function(x){sum(x^2)})
    if(any(is.infinite(prior*exp(tmp2))) || all(prior*exp(tmp2) == 0)){
      ind_min <- which.min(tmp2)
      result <- prior/prior[ind_min] * exp(tmp2 - tmp2[ind_min])
    }else{
      result <- prior * exp(tmp2)
    }
    result <- result/sum(result)
  })
  posterior <- t(posterior)
  class <- apply(posterior, 1, which.max)
  list(class = class, posterior = posterior)
}

# Variable screening.
var_screen <- function(x, y, num = NULL){
  n <- dim(x)[1]
  p <- dim(x)[2]
  count <- table(y)
  class <- names(count)
  K <- length(class)
  if(is.null(num)) num <- n
  screen <- function(x, y){
    mu_all <- mean(x)
    mu <- sapply(class, function(i){mean(x[y == i])})
    a <- sum((mu - mu_all)^2 * count)/(K - 1)
    tmp <- c()
    for (i in seq_len(K)) {
      tmp <- c(tmp, rep(mu[i], count[i]))
    }
    b_tmp <-  x - tmp
    b <- sum(b_tmp^2)/(n - K)
    a/b
  }
  f <- apply(x, 2, screen, y)
  ind <- order(f, decreasing = TRUE)[1:num]
  ind
}
