# ------------------------------
library(flexclust)
# Generate synthetic data with sparse and low-rank structure
Generate.data <- function(n, p, r, pi0, piK, K, sparsity = 0.02, mue = 1.5, sparsity.latent = 0.7, seed) {
  set.seed(seed)
  p.all <- p + r
  
  SS <- matrix(runif(n = p**2, min = -1) * rbinom(n = p**2, size = 1, prob = sparsity/2), ncol = p)
  SS[upper.tri(SS)] <- t(SS)[upper.tri(SS)]
  
  LL <- diag(r)
  SLX <- matrix(runif(n = p*r, min = -1) * rbinom(n = p*r, size = 1, prob = sparsity.latent), ncol = r, nrow = p)
  SS <- 0.5 * (SS + t(SS))
  
  Theta.tilde <- rbind(cbind(SS, SLX), cbind(t(SLX), LL))
  Theta.tilde <- Theta.tilde - min(eigen(Theta.tilde)$val) * diag(p.all) + 0.1 * diag(p.all)
  Theta.tilde <- cov2cor(Theta.tilde)
  
  L <- Theta.tilde[1:p, (p+1):p.all] %*% solve(Theta.tilde[(p+1):p.all, (p+1):p.all]) %*% Theta.tilde[(p+1):p.all, 1:p]
  S <- Theta.tilde[1:p, 1:p]
  Theta <- S - L
  
  # Generate data
  p.tidle=p/10
  Mu <- matrix(0, K, p)
  Mu[1, ] <- c(rep(-mue,p.tidle), rep(0, p-p.tidle))
  Mu[2, ] <- c(rep(-mue, p.tidle/2), rep(mue, p.tidle/2), rep(0, p-p.tidle))
  Mu[3, ] <- c(rep(mue, p.tidle), rep(0, p-p.tidle))
  
  data <- matrix(0, n, p)
  
  if (pi0 != 0) {
    member <- sample(0:K, n, replace = T, prob = c(pi0, piK))
    n.noise <- sum(member == 0)
    X <- cbind(matrix(runif(n.noise * (p.tidle/2), -20, 5), ncol =p.tidle/2), matrix(runif(n.noise * (p.tidle/2),-10,10), ncol =p.tidle/2))
    data[member == 0, ] <- cbind(X, mvrnorm(n.noise, rep(0, p - p.tidle), Theta[-(1:p.tidle), -(1:p.tidle)]))
  } else {
    member <- sample(1:K, n, replace = T, prob = piK)
  }
  
  for (k in 1:K) {
    data[member == k, ] <- mvrnorm(length(which(member == k)), Mu[k, ], solve(Theta))
  }
  
  return(list(
    data = data, S = S, L = L, Theta = Theta, 
    Mu = Mu, member = member
  ))
}


# ------------------------------
# Evaluate model performance using multiple metrics
Metric= function(data.test,memb.test,lambda1,lambda2,fit,S.true,Theta.true,S.hat,Theta.hat,L.hat,L.true,memb.true,memb.hat){
  p=ncol(S.true)
  FL.S = norm(S.hat- S.true,type="F")/norm(S.true,type="F")
  FL.L = norm(L.hat- L.true,type="F")/(norm(L.true,type="F")+0.001)
  FL.Theta = norm(Theta.hat- Theta.true,type="F")/norm(Theta.true,type="F")
  
  
  TPR =(sum((S.hat!=0) + (S.true!=0) == 2)-p)/(sum(S.true!=0)-p)
  FPR = sum((S.hat!=0) + (S.true==0) == 2)/sum(S.true==0)
  
  TPR.noise = sum((memb.hat==0) + (memb.true==0) == 2)/sum(memb.true==0)
  FPR.noise = sum((memb.hat==0) + (memb.true!=0) == 2)/sum(memb.true!=0)
  
  ind=which(memb.true!=0)
  ARI.c = comPart(memb.hat[ind],memb.true[ind], type="ARI")
  ARI.o=comPart(as.numeric(memb.hat!=0),as.numeric(memb.true!=0), type="ARI")
  cluster.pre=Predict(data.test,fit$pi,fit$mean,fit$Theta,G=3,icd=exp(fit$logicd),equalTheta=T)
  ARI.pre=comPart(cluster.pre,memb.test, type="ARI")
  
  
  rank=qr(L.hat)$rank
  index = as.data.frame(t(c(lambda1,lambda2,FL.S,FL.L,FL.Theta,rank,TPR, FPR,ARI.c,ARI.o,ARI.pre)))
  names(index) = c("lambda1","lambda2","FL.S","FL.L","FL.Theta","rank","TPR","FPR","ARI.c","ARI.o","ARI.pre")
  return(index)
}


# ------------------------------
# Predict cluster membership for test data
Predict <- function(data, pr, M, Theta, G, icd, equalTheta = TRUE) {
  N <- dim(data)[1L]
  P <- dim(data)[2L]
  
  gausscost <- (2 * pi)^(-P/2)
  
  psi <- rep(0, times = N)
  phi <- matrix(0, nrow = N, ncol = G + 1, byrow = TRUE, dimnames = NULL)
  
  for (j in 1:G) {
    Delta <- data - matrix(M[, j], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
    
    if (equalTheta) {
      smd <- .rowSums((Delta %*% Theta) * Delta, m = N, n = P, na.rm = TRUE)
      theta.det <- det(Theta)
    } else {
      smd <- .rowSums((Delta %*% Theta[,, j]) * Delta, m = N, n = P, na.rm = TRUE)
      theta.det <- det(Theta[,, j])
    }
    
    phi[, 1 + j] <- pr[1 + j] * gausscost * theta.det^(1/2) * exp(-0.5 * smd)
  }
  
  phi[, 1] <- pr[1] * icd
  psi <- .rowSums(phi, m = N, n = G + 1, na.rm = TRUE)
  tau <- phi * matrix(1/psi, nrow = N, ncol = G + 1, byrow = FALSE, dimnames = NULL)
  tau[is.na(tau)] <- 0
  
  cluster <- apply(tau, 1, which.max) - 1
  return(cluster)
}