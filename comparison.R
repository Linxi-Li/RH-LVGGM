# ==============================
library(glasso)
# (b) H-LVGGM: Heterogeneous LVGGM 
# Incorporates cluster-wise heterogeneity but ignores outlier noise components 
# lambda2.seq=1000 can force rank=0 to exclude latent variables
# ==============================
H.lvggm <- function(data,tau.ini,lambda.1,lambda.2,npr.max, erc, iter.max, tol) {
  N <- dim(data)[1L]
  P <- dim(data)[2L]
  G <- dim(tau.ini)[2L]
  
  tau.old <- tau.ini
  sumtau.old <- .colSums(tau.old, m = N, n = G, na.rm = TRUE)
  
  gausscost <- (2 * pi)^(-P/2)
  
  pr <- sumtau.new <- rep(0, times = G)
  psi <- tau.new <- wts <- rep(0, times = N)
  phi <- matrix(0, nrow = N, ncol = G, byrow = TRUE, dimnames = NULL)
  Delta <- data
  M.old <- M <- matrix(0, nrow = P, ncol = G, byrow = TRUE, dimnames = NULL)
  Theta.old <- Theta <- diag(1, P, P)
  dif <- tol + 1
  iter <- 0
  
  while ((iter <= iter.max) & (dif > tol)) {
    pr <- sumtau.old / N
    Xpp <- matrix(0, nrow = P, ncol = P, byrow = TRUE, dimnames = NULL)
    
    for (k in 1:G) {
      # Update mean vectors                       
      wts <- tau.old[, k] / sumtau.old[k]
      M[, k] <- .colSums(wts * data, m = N, n = P, na.rm = TRUE)
      Delta <- data - matrix(M[, k], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
      Xpp <- Xpp + crossprod(sqrt(wts) * Delta, y = NULL) * sumtau.old[k]
    }
    
    Xpp <- Xpp / sum(sumtau.old)
    E <- eigen(x = Xpp, symmetric = TRUE)
    D <- as.matrix(E$values)
    V <- E$vectors
    
    Dmin <- min(D)
    Dmax <- max(D)
    EigenRatio <- Dmax / Dmin
    
    if (EigenRatio > erc | !is.finite(EigenRatio) | any(D <= 0) | any(!is.finite(D))) {
      D <- GssERC(values = D, erc = erc, Dmin = Dmin, Dmax = Dmax, sumtau = sumtau.old[-1], P = P, G = G)
      Xpp <- V %*% (D[, 1] * t(V))
    }
    
    fit.Theta <- lrpsadmm(Xpp, lambda.1, lambda.2)
    Theta <- fit.Theta$S - fit.Theta$L
    
    # Update tau
    for (j in 1:G) {
      Delta <- data - matrix(M[, j], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
      smd <- .rowSums((Delta %*% Theta) * Delta, m = N, n = P, na.rm = TRUE)
      phi[, j] <- pr[j] * gausscost * (det(Theta))^(1/2) * exp(-0.5 * smd)
    }
    
    psi <- .rowSums(phi, m = N, n = G, na.rm = TRUE)
    tau.new <- phi * matrix(1/psi, nrow = N, ncol = G, byrow = FALSE, dimnames = NULL)
    tau.new[is.na(tau.new)] <- 0
    sumtau.new <- .colSums(tau.new, m = N, n = G, na.rm = TRUE)
    
    dif <- norm(M - M.old, type = "F") / norm(M.old, type = "F") + 
      norm(Theta - Theta.old, type = "F") / norm(Theta.old, type = "F")
    M.old <- M
    Theta.old <- Theta
    tau.old <- tau.new
    sumtau.old <- sumtau.new
    iter <- iter + 1
  }
  ans <- list(
    iter = iter,
    logicd = -Inf,
    pi = c(0, pr),
    mean = M,
    Theta = Theta,
    tau = tau.old,
    smd = smd,
    S = fit.Theta$S,
    L = fit.Theta$L,
    cluster = apply(tau.old, 1, rwhich.max),
    iloglik = sum(log(psi))
  )
  
  return(ans)
}


# ==============================
# Hyperparameter optimization for H-LVGGM
# ==============================
Compare.H.lvggm <- function(data, G, tau.ini, lambda1.seq, lambda2.seq, npr.max=0.5, erc=50, iter.max=100, tol=1e-03) {
  p <- ncol(data)
  n <- nrow(data)
  n.lambda1 <- length(lambda1.seq)
  n.lambda2 <- length(lambda2.seq)
  
  lambda.1 <- median(lambda1.seq)
  lambda.2 <- median(lambda2.seq)
  
  aAIC <- c()
  df <- c()
  
  for (l in 1:n.lambda2) {
    lambda.2 <- lambda2.seq[l]
    fit <- H.lvggm(data, tau.ini, lambda.1, lambda.2, npr.max, erc, iter.max, tol)
    df <- length(which(fit$S != 0)) + p * qr(fit$L)$rank
    aAIC[l] <- -2 * fit$iloglik + 2 * df
  }
  
  lambda.2 <- lambda2.seq[which.min(aAIC)]
  
  aAIC <- c()
  df <- c()
  result <- list()
  
  for (l in 1:n.lambda1) {
    lambda.1 <- lambda1.seq[l]
    result[[l]] <- fit <- H.lvggm(data, tau.ini, lambda.1, lambda.2, npr.max, erc, iter.max, tol)
    df <- length(which(fit$S != 0)) + p * qr(fit$L)$rank
    aAIC[l] <- -2 * fit$iloglik + 2 * df
  }
  
  lambda.1 <- lambda1.seq[which.min(aAIC)]
  result.final <- result[[which.min(aAIC)]]
  result.final$lambda.1 <- lambda.1
  result.final$lambda.2 <- lambda.2
  
  return(result.final)
}


# ==============================
# (d) Glasso: Classic Graphical Lasso 
# ==============================
Compare.glasso <- function(data, lambda1.seq) {
  p <- ncol(data)
  n <- nrow(data)
  n.lambda1 <- length(lambda1.seq)
  sigma <- cov(data)
  aAIC <- c()
  df <- c()
  result <- list()
  
  for (l in 1:n.lambda1) {
    lambda.1 <- lambda1.seq[l]
    result[[l]] <- fit <- glasso(s = sigma, rho = lambda.1)
    loglik <- (log(det(fit$wi)) - sum(diag(sigma %*% fit$wi))) * n
    df <- length(which(fit$wi != 0))
    aAIC[l] <- -2 * loglik + 2 * df
  }
  
  lambda.1 <- lambda1.seq[which.min(aAIC)]
  fit.final <- result[[which.min(aAIC)]]
  
  result.final <- list(
    Theta = fit.final$wi,
    S = fit.final$wi,
    L = matrix(0, p, p),
    cluster = rep(1, n),
    lambda.1 = lambda.1,
    mean = matrix(colMeans(data), p)
  )
  
  return(result.final)
}


# ==============================
# (e) OR-Glasso: Oracle Glasso with True Labels 
# ==============================
Compare.or.glasso <- function(data, K, member, lambda1.seq) {
  p <- ncol(data)
  n <- nrow(data)
  Theta <- array(0, dim = c(p, p, K))
  Mu <- matrix(0, p, K)
  pi <- rep(0, K)
  
  for (k in 1:K) {
    data.k <- data[member == k, ]
    fit <- Compare.glasso(data.k, lambda1.seq)
    Theta[,, k] <- fit$Theta
    Mu[, k] <- colMeans(data[member == k, ])
    pi[k] <- sum(member == k) / n
  }
  
  Theta[abs(Theta) < 1e-3] <- 0
  
  result <- list(
    K = K,
    mean = Mu,
    Theta = Theta,
    cluster = member,
    pi = pi,
    S = Theta,
    L = matrix(0, p, p)
  )
  
  return(result)
}


# ==============================
# (c) OR-LVGGM: Oracle LVGGM with True Labels 
# ==============================
Compare.or.lvggm <- function(data, K, member, lambda1.seq, lambda2.seq) {
  p <- ncol(data)
  n <- nrow(data)
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  L <- array(0, dim = c(p, p, K))
  Mu <- matrix(0, p, K)
  pi <- rep(0, K)
  
  for (k in 1:K) {
    data.k <- data[member == k, ]
    fit <- Compare.lvggm(data.k, lambda1.seq, lambda2.seq)
    Theta[,, k] <- fit$Theta
    S[,, k] <- fit$S
    L[,, k] <- fit$L
    Mu[, k] <- colMeans(data[member == k, ])
    pi[k] <- sum(member == k) / n
  }
  
  Theta[abs(Theta) < 1e-3] <- 0
  
  result <- list(
    K = K,
    mean = Mu,
    S = S,
    L = L,
    Theta = Theta,
    cluster = member,
    pi = pi
  )
  
  return(result)
}


# ==============================
# (a) LVGGM: Classic Latent Variable Graphical Model 
# Homogeneous model without outlier consideration
# ==============================
Compare.lvggm <- function(data, lambda1.seq, lambda2.seq) {
  p <- ncol(data)
  n <- nrow(data)
  sigma <- cov(data)
  n.lambda1 <- length(lambda1.seq)
  n.lambda2 <- length(lambda2.seq)
  
  lambda.1 <- median(lambda1.seq)
  lambda.2 <- median(lambda2.seq)
  
  aAIC <- c()
  df <- c()
  
  for (l in 1:n.lambda2) {
    lambda.2 <- lambda2.seq[l]
    fit <- lrpsadmm(sigma, lambda.1, lambda.2)
    Theta <- fit$S - fit$L
    loglik <- (log(det(Theta) + 1e-6) - sum(diag(sigma %*% Theta))) * n
    df <- length(which(fit$S != 0)) + p * qr(fit$L)$rank
    aAIC[l] <- -2 * loglik + 2 * df
  }
  
  lambda.2 <- lambda2.seq[which.min(aAIC)]
  
  aAIC <- c()
  df <- c()
  result <- list()
  
  for (l in 1:n.lambda1) {
    lambda.1 <- lambda1.seq[l]
    result[[l]] <- fit <- lrpsadmm(sigma, lambda.1, lambda.2)
    Theta <- fit$S - fit$L
    loglik <- (log(det(Theta) + 1e-6) - sum(diag(sigma %*% Theta))) * n
    df <- length(which(fit$S != 0)) + p * qr(fit$L)$rank
    aAIC[l] <- -2 * loglik + 2 * df
  }
  
  lambda.1 <- lambda1.seq[which.min(aAIC)]
  result.final <- result[[which.min(aAIC)]]
  result.final$lambda.1 <- lambda.1
  result.final$lambda.2 <- lambda.2
  result.final$Theta <- result.final$S - result.final$L
  result.final$cluster <- rep(1, n)
  result.final$mean <- colMeans(data)
  
  return(result.final)
}


# ==============================
# Cluster initialization without noise consideration
# ==============================
InitClust <- function(data, G) {
  hc <- hclust(dist(data, method = "euclidean"), method = "ward.D")
  cluster <- cutree(hc, k = G)
  cl <- sort(unique(cluster))
  G <- length(cl)
  N <- length(cluster)
  A <- matrix(0, nrow = N, ncol = G)
  
  for (j in 1:G) {
    A[cluster == cl[j], j] <- 1
  }
  
  return(A)
}