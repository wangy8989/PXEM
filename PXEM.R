#input X, y, Z
PXEM <- function(X, y, Z, max.iter = 3000 , eps = 0.000001){
  delta <- 1L
  # X normalization
  
  p <- dim(X)[[2]]
  n <- dim(X)[[1]]
  X <- (X - colMeans(X))/(resample::colStdevs(X))*sqrt(p)
  
  
  # initialzation 
  iter <- 0
  
  # ? population variance
  w <- MASS::ginv((t(Z) %*% Z)) %*% t(Z) %*% y
  Sig_e2 <- var(y-Z %*% w)[1,1]
  Sig_b2 <- var(y-Z %*% w)[1,1]
  
  while(iter < max.iter | abs(elbo) > eps){
    # E-step
    MSig_inv <- (1/Sig_e2) * t(X) %*% X + (1/Sig_b2) * diag(p)
    MSig <- MASS::ginv((1/Sig_e2) * t(X) %*% X + (1/Sig_b2) * diag(p))
    mu <- MSig %*% ((1/Sig_e2) * t(X) %*% (y-Z %*% w))
    
    elbo <- - n/2 * log(2 * pi * Sig_e2[1,1]) - p/2 * log(2 * pi * Sig_b2[1,1]) - 0.5 * (1/Sig_e2[1,1]) * t(y - Z %*% w- delta * X %*% mu) %*%
      (y - Z %*% w- delta * X %*% mu) - 0.5 * (1/Sig_b2[1,1]) * t(mu) %*% mu -sum(diag(( delta * delta * 0.5 * (1/Sig_e2[1,1]) * t(X) %*% X + 
                                                                                  0.5 * (1/Sig_b2[1,1]) * diag(p) ) %*% MSig))
      
    delta <- (t(y-Z %*% w) %*% X %*% mu)/(t(mu) %*% t(X) %*% X %*% mu+sum(diag(t(X) %*% X %*% MSig)))
    w <- MASS::ginv(t(Z) %*% Z) %*% t(Z) %*% (y-delta * X %*% mu)
    Sig_e2 <- (t(y-delta * X %*% mu-Z %*% w) %*% (y-delta * X %*% mu-Z %*% w)+delta*delta*sum(diag(t(X) %*% X %*% MSig)))/n
    Sig_b2 <- (t(mu) %*% mu + sum(diag(MSig)))/p
    
    # reset delta=1
    delta <- 1L
    iter <- iter + 1
  }
 
}

