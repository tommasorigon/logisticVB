logit_CAVI <- function(X, y, prior, tol = 1e-16, maxiter=10000){
  
  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")
  
  # Compute the log-determinant of a matrix
  ldet <- function(X) {
    if(!is.matrix(X)) return(log(X))
    determinant(X,logarithm = TRUE)$modulus
  }

  lowerbound <- numeric(maxiter)
  p          <- ncol(X)
  
  P    <- solve(prior$Sigma)
  mu   <- prior$mu
  Pmu  <- c(P%*%mu)
  Pdet <- ldet(P)
  
  # Initialization for omega equal to 0.25
  P_vb       <- crossprod(X*rep(1/4,n),X) + P
  Sigma_vb   <- solve(P_vb) 
  mu_vb      <- Sigma_vb %*% (crossprod(X,y - 0.5) + Pmu)
  eta        <- c(X%*%mu_vb)
  xi         <- sqrt(eta^2 +  rowSums(X %*% Sigma_vb * X))
  omega      <- tanh(xi/2)/(2*xi); 
  omega[is.nan(omega)] <- 0.25
  
  lowerbound[1]  <- -0.5*p + 0.5*ldet(Sigma_vb) + 0.5*Pdet - 0.5*t(mu_vb - mu)%*%P%*%(mu_vb - mu) + sum((y-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(P %*% Sigma_vb))
  
  # Iterative procedure
  for(t in 2:maxiter){
  
    P_vb       <- crossprod(X*omega,X) + P
    Sigma_vb   <- solve(P_vb) 
    mu_vb      <- Sigma_vb %*% (crossprod(X,y-0.5) + Pmu)
    
    #Update of xi
    eta        <- c(X%*%mu_vb)
    xi         <- sqrt(eta^2 +  rowSums(X %*% Sigma_vb * X))
    omega      <- tanh(xi/2)/(2*xi); 
    omega[is.nan(omega)] <- 0.25
    
    lowerbound[t]  <- -0.5*p + 0.5*ldet(Sigma_vb) + 0.5*Pdet - 0.5*t(mu_vb - mu)%*%P%*%(mu_vb - mu) + sum((y-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(P %*% Sigma_vb))
    
    if(abs(lowerbound[t] - lowerbound[t-1]) < tol) return(list(mu = matrix(mu_vb,p,1), Sigma=matrix(Sigma_vb,p,p), 
                                                               Convergence=cbind(Iteration=(1:t)-1, 
                                                                                 Lowerbound=lowerbound[1:t]), xi=xi))
  }
  stop("The algorithm has not reached convergence")
}


logit_SVI <- function(X, y, prior, iter, tau, kappa){
  
  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")
  
  p          <- ncol(X)

  P    <- solve(prior$Sigma)
  mu   <- prior$mu
  Pmu  <- c(P%*%mu)

  # Initialization of Eta1 and Eta2
  Eta1       <- Pmu
  Eta2       <- P

  Eta1_out <- Eta1
  Eta2_out <- Eta2
  
  # Iterative procedure
  for(t in 1:iter){
  
    # Sample the observation
    id  <- sample.int(n,1)
    x_i <- X[id,]
    y_i <- y[id]
    
    # Update the local parameter
    Sigma_vb   <- solve(Eta2_out) 
    mu_vb      <- Sigma_vb%*%Eta1_out
    
    eta_i   <- c(crossprod(x_i, mu_vb))
    xi_i    <- sqrt(eta_i^2 +  rowSums(x_i %*% Sigma_vb * x_i))
    omega_i <- tanh(xi_i/2)/(2*xi_i); 
    
    Eta1       <- n*x_i*(y_i-0.5) + Pmu
    Eta2       <- n*tcrossprod(x_i*omega_i,x_i) + P
    
    # Update the final estimates
    rho      <- 1/(t + tau)^kappa 
    Eta1_out <- (1 - rho)*Eta1_out + rho*Eta1
    Eta2_out <- (1 - rho)*Eta2_out + rho*Eta2
  }
  
  # Output
  Sigma_vb   <- matrix(solve(Eta2_out),p,p)
  mu_vb      <- matrix(Sigma_vb%*%Eta1_out,p,1)
  
  return(list(mu=mu_vb,Sigma=Sigma_vb))
}