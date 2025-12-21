# ------------------------------
library(MASS)
# Update A matrix with positive definiteness projection
UpdateA <- function(Sigma, S, L, U, mu, epsilon = 0.001) {
  X1 <- mu * (S - L) - Sigma - U
  X2 <- X1 %*% X1 + 4 * mu * diag(dim(X1)[1])
  eig <- eigen(X2, symmetric = TRUE)
  sqrtX2 <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  A <- (X1 + sqrtX2) / (2 * mu)
  A <- (A + t(A)) / 2
  eig.A <- eigen(A, symmetric = TRUE)
  new.eigenvalues <- pmax(eig.A$values, epsilon)
  A.proj <- eig.A$vectors %*% diag(new.eigenvalues) %*% t(eig.A$vectors)
  A.proj<-apply(A.proj, 1:2, as.numeric)
  return(A.proj)
}

# ------------------------------
# Update L matrix with non-negativity projection
UpdateL <- function(A, S, U, l2, mu) {
  lp2 <- l2 / mu
  X1 <- S - A - (U / mu)
  eig <- eigen(X1, symmetric = TRUE) 
  eigVal <- MCP.op(eig$values, lp2)
  L.temp <- eig$vectors %*% diag(eigVal) %*% t(eig$vectors)
  L.temp <- (L.temp + t(L.temp)) / 2
  eig.L <- eigen(L.temp, symmetric = TRUE)
  new.eigenvalues <- pmax(eig.L$values, 0) 
  L.proj <- eig.L$vectors %*% diag(new.eigenvalues) %*% t(eig.L$vectors)
  L.proj<-apply(L.proj, 1:2, as.numeric)
  return(L.proj)
}

# ------------------------------
# MCP (Minimax Concave Penalty) operator
MCP.op = function(z,lambda,a=3){
  return( Soft.op(z,lambda)/(1-1/a) * (abs(z) - a*lambda <= 0) + z * (abs(z) - a*lambda > 0) )
}

# ------------------------------
# Soft thresholding operator
Soft.op = function(z,lambda){
  return((abs(z) - lambda)*(abs(z) - lambda > 0)*ifelse(Mod(z) == 0, 0, z / Mod(z)))
}

# ------------------------------
UpdateS <- function(A, L, U, l1, mu) {
  lp1 <- l1 / mu
  X1 <- A + L + (U / mu)
  S <- MCP.op(X1,lp1)
  S <-apply(S, 1:2, as.numeric)
  return(S)
}

# ------------------------------
UpdateU <- function(A, S, L, U, mu) {
  U <- U + mu * (A - S + L)
  return(U)
}

# ------------------------------
# ADMM-based optimization for LVGGM
lrpsadmm <- function(Sigma, lambda.1, lambda.2,init=NULL,maxiter=100, mu=1, abs.tol=1e-03, rel.tol=1e-02) {
  
  p <- dim(Sigma)[1]
  
  if (is.null(init)) {
    S <- diag(p)
    L <- S * 0.0
    U <- S * 0.0
  } else {
    parameters <- init
    S <- init$S
    L <- init$L
    U <- init$U 
  }
  
  S <- S 
  history <- matrix(NA, nrow=0, ncol=3)
  colnames(history) <- c('Iteration', 'Objval','r.norm')
  parameters <- list()
  parameters$termcode <- -1
  parameters$termmsg <- "Maximum number of iterations reached."
  
  for (i in 1:maxiter) {
    
    # Update A
    A <- UpdateA(Sigma, S, L, U, mu)
    
    # Update S
    S.old <- S
    S <- UpdateS(A, L, U, lambda.1, mu)
    
    if(any(diag(S) == 0)) {
      S <- diag(p)
      L <- S * 0
      U <- S * 0
      parameters$termcode <- -2
      parameters$termmsg <- "Shrinkage too strong: sparse component is empty."
      break()
    }
    
    # Update L
    L.old <- L
    L <- UpdateL(A, S, U, lambda.2, mu)
    
    # Update U
    U <- UpdateU(A, S, L, U, mu)
    
    # Diagnostics
    objval <- obj.func(Sigma, A, S, L, lambda.1, lambda.2)
    r.norm <- norm(A - (S-L), 'F')
    history <- rbind(history, c(i,objval,r.norm))
    
    if(  (r.norm < abs.tol) ) {
      break()
    }
  }
  
  parameters$S <- S
  parameters$L <- L
  parameters$U <- U
  parameters$history <- as.data.frame(history)
  return(parameters)
}

# ------------------------------
# Objective function for LVGGM
obj.func <- function(Sigma, A, S, L, l1, l2) {
  evals <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  evals[abs(evals) < 1e-09] <- 1e-8
  
  if (any(evals < 0)) {
    return(NaN)
  }
  
  ll <- sum(diag(Sigma %*% A)) - sum(log(evals))
  ll <- ll + l1 * sum(abs(S)) + l2 * sum(diag(L))
  return(ll)
}

# ------------------------------
GssERC <- function(values, erc, Dmin, Dmax, sumtau, P, G) {
  tol <- sqrt(.Machine$double.eps)
  iter.max <- 99
  iter <- 0L
  grc <- {
    sqrt(5) - 1
  }/2
  a <- log(max(Dmin, .Machine$double.eps))
  b <- log(min(Dmax, .Machine$double.xmax))
  tD <- values
  l <- a
  x1 <- a + {
    1 - grc
  } * {
    b - a
  }
  x2 <- a + grc * {
    b - a
  }
  
  iter <- iter + 1L
  l <- exp(x1)
  tD <- values
  tD[values < l] <- l
  tD[values > erc * l] <- l * erc
  f1 <- sum(sumtau * .colSums(log(tD) + values/tD, m = P, n = 1, na.rm = FALSE))
  
  iter <- iter + 1L
  l <- exp(x2)
  tD <- values
  tD[values < l] <- l
  tD[values > erc * l] <- l * erc
  f2 <- sum(sumtau * .colSums(log(tD) + values/tD, m = P, n = 1, na.rm = FALSE))
  
  while ({
    abs(b - a) > tol
  } & {
    iter < iter.max
  }) {
    if (f1 > f2) {
      a <- x1
      x1 <- x2
      x2 <- a + grc * {
        b - a
      }
      f1 <- f2
      
      iter <- iter + 1L
      l <- exp(x2)
      tD <- values
      tD[values < l] <- l
      tD[values > erc * l] <- l * erc
      f2 <- sum(sumtau * .colSums(log(tD) + values/tD, m = P, n = 1, na.rm = FALSE))
    } else {
      b <- x2
      x2 <- x1
      x1 <- a + {
        1 - grc
      } * {
        b - a
      }
      f2 <- f1
      
      l <- exp(x1)
      iter <- iter + 1L
      tD <- values
      tD[values < l] <- l
      tD[values > erc * l] <- l * erc
      f1 <- sum(sumtau * .colSums(log(tD) + values/tD, m = P, n = 1, na.rm = FALSE))
    }
  }
  
  l <- exp({
    a + b
  }/2)
  tD <- values
  tD[values < l] <- l
  tD[values > erc * l] <- l * erc
  return(tD)
}

# ------------------------------
EqNPC <- function(x, icd, sumtau.old, pr, psi, N, npr.max) {
  if (icd * x == 0) {
    return(-N * npr.max)
  } else {
    return(sum(icd * x/{
      icd * x + {
        {
          1 - x
        }/{
          N - sumtau.old[1]
        }
      } * {
        {
          psi - pr[1] * icd
        } * N
      }
    }) - N * npr.max)
  }
}

# ------------------------------
rwhich.min <- function(x) {
  xmin <- min(x, na.rm = TRUE)
  tmp <- which(x == xmin)
  
  if (length(tmp) == 1) {
    return(tmp)
  } else {
    return(sample(tmp, size = 1, replace = FALSE))
  }
}

# ------------------------------
CountUniqueRows <- function(x) {
  if (is.vector(x)) {
    if (length(x) == 0) {
      return(0)
    } else {
      return(length(unique(x)))
    }
  } else {
    n <- nrow(x)
    if (n == 0) {
      return(0)
    } else {
      if (ncol(x) == 1) {
        return(length(unique(x[, 1])))
      } else {
        sx <- x[do.call(order, as.list(as.data.frame(x))), ]
        flags <- rep(0L, times = n)
        
        for (i in 2:n) {
          flags[i] <- prod(sx[i, ] == sx[i - 1, ])
        }
        
        return({
          n - sum(flags)
        })
      }
    }
  }
}

# ------------------------------
Cluster2Assign <- function(cluster) {
  cl <- sort(unique(cluster[cluster != 0]))
  G <- length(cl)
  N <- length(cluster)
  A <- matrix(0, nrow = N, ncol = 1 + G)
  
  if (sum(cluster == 0) > 0) {
    A[cluster == 0, 1] <- 1
  }
  
  for (j in 1:G) {
    A[cluster == cl[j], j + 1] <- 1
  }
  
  return(A)
}

# ------------------------------
rwhich.max <- function(x) {
  xmax <- max(x, na.rm = TRUE)
  tmp <- which(x == xmax)
  
  if (length(tmp) == 1) {
    return(tmp)
  } else {
    return(sample(tmp, size = 1, replace = FALSE))
  }
}

# ------------------------------
wecdf <- function(x, weights) {
  sw <- sum(weights)
  
  if (sw == 0) {
    ans <- rep(0, length(x))
  } else {
    ans <- wecdf.aux(x, weights/sw)(x)
  }
  
  return(ans)
}

# ------------------------------
wecdf.aux <- function(x, weights) {
  ox <- order(x)
  x <- x[ox]
  ow <- weights[ox]
  n <- length(x)
  
  rval <- approxfun(x, cumsum(ow), method = "constant", yleft = 0, yright = 1, 
                    f = 0, ties = "ordered")
  
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}

# ------------------------------
# Main model training using EM algorithm with robust and heterogeneous components
RH.lvggm<- function(data,tau.ini,lambda.1,lambda.2,npr.max, erc, iter.max, tol) {
  N <- dim(data)[1L]
  P <- dim(data)[2L]
  G <- dim(tau.ini)[2L] - 1
  
  tau.old <- tau.ini
  sumtau.old <- .colSums(tau.old, m = N, n = G + 1, na.rm = TRUE)
  
  gausscost <- (2 * pi)^(-P/2)
  
  pr <- sumtau.new <- rep(0, times = 1 + G)
  psi <- tau.new <- wts <- rep(0, times = N)
  phi <- matrix(0, nrow = N, ncol =  G + 1, byrow = TRUE, dimnames = NULL)
  Delta <- data
  M.old=M <- matrix(0, nrow = P, ncol = G, byrow = TRUE, dimnames = NULL)
  Theta.old=Theta=diag(1,P,P)
  em.failed <- FALSE
  flag <- rep(FALSE, times = 4)
  dif=tol+1
  iter=0
  while((iter<=iter.max)&(dif>tol)) {
    pr <- sumtau.old/N
    Xpp <- matrix(0, nrow = P, ncol = P, byrow = TRUE, dimnames = NULL)
    for (k in 1:G) {
      # update mean vectors                       
      wts <- tau.old[, 1 + k]/sumtau.old[1 + k]
      M[, k] <- .colSums(wts * data, m = N, n = P, na.rm = TRUE)
      Delta <-  data - matrix(M[, k], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
      Xpp <- Xpp+crossprod(sqrt(wts) * Delta, y = NULL)*sumtau.old[1 + k]
    }
    
    Xpp<-Xpp/sum(sumtau.old[-1])
    
    E <- eigen(x = Xpp, symmetric = TRUE)
    D <- as.matrix(E$values)
    V<- E$vectors
    
    Dmin <- min(D)
    Dmax <- max(D)
    EigenRatio <- Dmax/Dmin
    if ({
      EigenRatio > erc
    } | !is.finite(EigenRatio) | any(D <= 0) | any(!is.finite(D))) {
      D <- GssERC(values =D, erc = erc, Dmin = Dmin,Dmax = Dmax, sumtau = sumtau.old[-1],P = P, G = G)
      Xpp <- V%*% (D[,1]* t(V))
    }
    
    #update theta
    fit.Theta<- lrpsadmm(Xpp, lambda.1, lambda.2)
    
    Theta=fit.Theta$S-fit.Theta$L
    
    cc=c()
    for(j in 1:P){
      c2=sum(tau.old[,1]/pr[1]*data[,j]^2)/N
      c1=sum(tau.old[,1]/pr[1]*data[,j])/N
      cc[j]=2*sqrt(3*(c2-c1^2))
    }
    icd=1/prod(cc)
    #update tau
    for (j in 1:G) {
      Delta <- data - matrix(M[, j], nrow = N, ncol = P, byrow = TRUE, dimnames = NULL)
      smd <- .rowSums((Delta %*% Theta)* Delta, m = N, n = P, na.rm = TRUE)
      phi[, 1 + j] <- pr[1 + j] * gausscost *(det(Theta))^(1/2)  * (exp(-0.5 * smd))
    }
    phi[, 1] <- pr[1] * icd
    psi <- .rowSums(phi, m = N, n = G + 1, na.rm = TRUE)
    tau.new <- phi * matrix(1/psi, nrow = N, ncol = 1 + G, byrow = FALSE, dimnames = NULL)
    sumtau.new <- .colSums(tau.new, m = N, n = G + 1, na.rm = TRUE)
    
    
    if (any(!is.finite(tau.new)) | any(sumtau.new[-1] == 0)) {
      em.failed <- TRUE
      flag[1] <- TRUE
    }
    if (!em.failed) {
      tmp <- (.rowSums(tau.new, m = N, n = {
        G + 1
      }, na.rm = TRUE) == 0)
      if (any(tmp)) {
        tau.new[tmp, ] <- 0
        if (icd > 0) {
          tau.new[tmp, 1] <- 1
        }
        if (icd == 0) {
          i.tmp <- which(tmp)
          n.tmp <- length(i.tmp)
          ED <- matrix(0, nrow = n.tmp, ncol = G, byrow = TRUE, dimnames = NULL)
          for (j in 1:G) {
            tmpDelta <- data[i.tmp, ] - matrix(M[, j], byrow = TRUE, nrow = n.tmp, 
                                               ncol = P, dimnames = NULL)
            ED[, j] <- .rowSums(tmpDelta * tmpDelta, m = n.tmp, n = P, na.rm = TRUE)
          }
          cl.tmp <- 1 + apply(-ED, 1, rwhich.max)
          for (k in 1:n.tmp) {
            tau.new[i.tmp[k], cl.tmp[k]] <- 1
          }
        }
        sumtau.new <- .colSums(tau.new, m = N, n = {
          G + 1
        }, na.rm = TRUE)
      }
    }
    if (!em.failed & icd > 0 & {
      sumtau.new[1]/N > npr.max
    }) {
      z <- NULL
      try(z <- uniroot(f =EqNPC, icd = icd, sumtau.old = sumtau.old, pr = pr, 
                       psi = psi, N = N, npr.max = npr.max, lower = 0, upper = npr.max, 
                       tol = .Machine$double.eps, maxiter = 500, trace = 0), silent = TRUE)
      if (is.null(z)) {
        flag[2] <- TRUE
        em.failed <- TRUE
      }
      else {
        flag[3] <- TRUE
        for (j in 1:G) {
          phi[, 1 + j] <- phi[, 1 + j]/pr[1 + j]
        }
        pr <- c(z$root, {
          {
            1 - z$root
          }/{
            N - sumtau.old[1]
          }
        } * sumtau.old[-1])
        for (j in 1:G) {
          phi[, 1 + j] <- pr[1 + j] * phi[, 1 + j]
        }
        phi[, 1] <- pr[1] * icd
        psi <- .rowSums(phi, m = N, n = {
          G + 1
        }, na.rm = TRUE)
        tau.new <- phi * matrix(1/psi, nrow = N, ncol = {
          1 + G
        }, byrow = FALSE, dimnames = NULL)
        sumtau.new <- .colSums(tau.new, m = N, n = {
          G + 1
        }, na.rm = TRUE)
        anpr <- mean({
          apply(tau.new, 1, rwhich.max) - 1
        } == 0)
        if (anpr > {
          npr.max + 0.01
        }) {
          flag[2] <- TRUE
          em.failed <- TRUE
        }
        
      }
    }
    dif=norm(M-M.old,type = "F")/norm(M.old,type = "F")+norm(Theta-Theta.old,type = "F")/norm(Theta.old,type = "F")
    M.old <- M
    Theta.old <- Theta
    tau.old <- tau.new
    sumtau.old <- sumtau.new
    iter=iter+1
  }
  ans=list()
  ans$iter <- iter
  ans$logicd <- log(icd)
  ans$pi <- pr
  ans$mean <- M
  ans$Theta <- Theta
  ans$tau <- tau.old
  ans$smd <- smd
  ans$S<- fit.Theta$S
  ans$L<- fit.Theta$L
  ans$cluster <- apply(tau.old, 1, rwhich.max) - 1
  ans$iloglik<- sum(log(psi))
  return(ans)
}
#Hyperparameter optimization using Akaike Information Criterion (AIC)
OT.RH.lvggm<-function(data,G,tau.ini,lambda1.seq,lambda2.seq,npr.max=0.5, erc=20,iter.max=100, tol=1e-03)
{
  p=ncol(data)
  n=nrow(data)
  
  n.lambda1<- length(lambda1.seq)
  n.lambda2<- length(lambda2.seq)
  
  
  
  lambda.1 = median(lambda1.seq)
  lambda.2 = median(lambda2.seq)
  
  
  aAIC=c()
  df=c()
  for (l in 1:n.lambda2){
    lambda.2=lambda2.seq[l]
    fit=RH.lvggm(data, tau.ini,lambda.1,lambda.2,npr.max, erc,  iter.max, tol)
    df= length(which(fit$S!= 0))+ p*qr(fit$L)$rank
    aAIC[l] =  -2 * fit$iloglik+2*df
  }
  lambda.2 = lambda2.seq[which.min(aAIC)]
  
  aAIC=c()
  df=c()
  result=list()
  for (l in 1:n.lambda1){
    lambda.1=lambda1.seq[l]
    result[[l]] =fit=RH.lvggm(data, tau.ini,lambda.1,lambda.2,npr.max, erc,iter.max, tol)
    df= length(which(fit$S!= 0))+ p*qr(fit$L)$rank
    aAIC[l] =  -2 * fit$iloglik+2*df
  }
  lambda.1 = lambda1.seq[which.min(aAIC)]
  result.final=result[[which.min(aAIC)]]
  result.final$lambda.1=lambda.1
  result.final$lambda.2=lambda.2
  return(result.final)
}
# Initialize cluster assignments with noise point detection
InitClust.nosie<- function(data,G) {
  hc <- hclust(dist(data,method = "euclidean"),method = "ward.D")
  cluster<- cutree(hc,k=G+1)
  label.counts <- table(cluster)
  min.count <- min(label.counts)
  cluster[cluster == which(label.counts == min.count)] <- 0
  cl <- sort(unique(cluster[cluster != 0]))
  G <- length(cl)
  N <- length(cluster)
  A <- matrix(0, nrow = N, ncol = 1 + G)
  A[cluster == 0, 1] <- 1
  for (j in 1:G) {
    A[cluster == cl[j], j + 1] <- 1
  }
  return(A)
}