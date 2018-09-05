logit_VB <- function(x, y,  tol = 1e-16, prior, maxiter=10000){
  
  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  
  lowerbound <- numeric(maxiter)
  p          <- ncol(x)
  
  P    <- solve(prior$Sigma)
  mu   <- prior$mu
  Pmu  <- c(P%*%mu)
  Pdet <- ldet(P)
  
  # Initialization for omega equal to 0.25
  P_vb       <- crossprod(x*rep(1/4,n),x) + P
  Sigma_vb   <- solve(P_vb) 
  mu_vb      <- Sigma_vb %*% (crossprod(x,y - 0.5) + Pmu)
  eta        <- c(x%*%mu_vb)
  xi         <- sqrt(eta^2 +  rowSums(x %*% Sigma_vb * x))
  omega      <- tanh(xi/2)/(2*xi); 
  omega[is.nan(omega)] <- 0.25
  
  lowerbound[1]  <- -0.5*p + 0.5*ldet(Sigma_vb) + 0.5*Pdet - 0.5*t(mu_vb - mu)%*%P%*%(mu_vb - mu) + sum((y-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(P %*% Sigma_vb))
  
  # Iterative procedure
  for(t in 2:maxiter){
  
    P_vb       <- crossprod(x*omega,x) + P
    Sigma_vb   <- solve(P_vb) 
    mu_vb      <- Sigma_vb %*% (crossprod(x,y-0.5) + Pmu)
    
    #Update of xi
    eta        <- c(x%*%mu_vb)
    xi         <- sqrt(eta^2 +  rowSums(x %*% Sigma_vb * x))
    omega      <- tanh(xi/2)/(2*xi); 
    omega[is.nan(omega)] <- 0.25
    
    lowerbound[t]  <- -0.5*p + 0.5*ldet(Sigma_vb) + 0.5*Pdet - 0.5*t(mu_vb - mu)%*%P%*%(mu_vb - mu) + sum((y-0.5)*eta +log(plogis(xi)) - 0.5*xi) - 0.5*sum(diag(P %*% Sigma_vb))
    
    if(abs(lowerbound[t] - lowerbound[t-1]) < tol) return(list(mu=mu_vb, Sigma=Sigma_vb, 
                                                               Convergence=cbind(Iteration=(1:t)-1, 
                                                                                 Lowerbound=lowerbound[1:t]), xi=xi))
  }
  stop("The algorithm has not reached convergence")
}


logit_SVB <- function(x, y,  tol = 1e-16, prior, maxiter=1000, tau=1, kappa=0.7){
  
  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  
  p          <- ncol(x)

  P    <- solve(prior$Sigma)
  mu   <- prior$mu
  Pmu  <- c(P%*%mu)

  # Initialization of Eta1 and Eta2
  Eta1       <- Pmu
  Eta2       <- P

  Eta1_out <- Eta1
  Eta2_out <- Eta2
  
  # Iterative procedure
  for(t in 1:maxiter){
  
    # Sample the observation
    id  <- sample.int(n,1)
    x_i <- x[id,]
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
  Sigma_vb   <- solve(Eta2_out) 
  mu_vb      <- Sigma_vb%*%Eta1_out
  
  return(list(mu=mu_vb,Sigma=Sigma_vb))
}


# For internal use only
ldet <- function(x) {
  if(!is.matrix(x)) return(log(x))
  determinant(x,logarithm = TRUE)$modulus
}