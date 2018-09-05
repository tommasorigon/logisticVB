#' IRLS algorithm for logistic regression
#'
#'@importFrom Rcpp evalCpp sourceCpp
#'@useDynLib logistic
#'
#'
#' @param x  A \verb{n x p} dimensional design matrix. The intercept must be manually included.
#' @param y  A \verb{n} dimensional binary vector.
#' @param tol Threshold at which the algorithm stops.
#' @param beta_start Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients.
#' @param maxiter Maximum number of iteration. If reached, the algorithm raises an error.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
#' 
#' @export
#' 
logit_NR <- function(x, y,  tol = 1e-16, beta_start = NULL, maxiter=100){

  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")

  loglik <- numeric(maxiter)

  # Initialization (If null, implicitely initialized at beta=0)
  if(is.null(beta_start)) {beta <- solve(crossprod(x/4,x),crossprod(x,y-0.5))} else {beta <- beta_start}
  # Initialization
  eta       <- c(x%*%beta)
  prob      <- 1/(1+exp(-eta))
  w         <- prob*(1-prob)

  # First value of the likelihood
  loglik[1] <- sum(y*eta - log(1+exp(eta)))

  # Iterative procedure
  for(t in 2:maxiter){
    beta       <- beta + solve(qr(crossprod(x*w,x)),crossprod(x,y-prob))
    eta        <- c(x%*%beta)
    prob       <- 1/(1+exp(-eta))
    w          <- prob*(1-prob)
    loglik[t]  <- sum(y*eta - log(1+exp(eta)))
    if(loglik[t] - loglik[t-1] < tol) return(list(beta=beta,vcov=solve(crossprod(x*w,x)),
                                                  Convergence=cbind(Iteration=(1:t)-1,
                                                                    Loglikelihood=loglik[1:t])))
  }
  stop("The algorithm has not reached convergence")
}

#' EM algorithm for logistic regression
#'
#'
#' @param x  A \verb{n x p} dimensional design matrix. The intercept must be manually included.
#' @param y  A \verb{n} dimensional binary vector.
#' @param tol Threshold at which the algorithm stops.
#' @param beta_start Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients.
#' @param maxiter Maximum number of iteration. If reached, the algorithm raises an error.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
#' 
#' @export
#' 
logit_EM <- function(x, y, tol = 1e-16, beta_start = NULL, maxiter=10000){

  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")

  loglik <- numeric(maxiter)
  Xy     <- crossprod(x,y-0.5)

  # Initialization (If null, implicitely initialized at beta=0)
  if(is.null(beta_start)) {beta <- solve(crossprod(x/4,x),crossprod(x,y-0.5))}else{beta <- beta_start}

  # Initialization
  eta        <- c(x%*%beta)
  w          <- c(EM_weights(eta))


  # First value of the likelihood
  loglik[1] <- sum(y*eta - log(1+exp(eta)))

  # Iterative procedure
  for(t in 2:maxiter){
    beta       <- solve(qr(crossprod(x*w,x)),Xy)
    eta        <- c(x%*%beta)
    w          <- c(EM_weights(eta))
    loglik[t]  <- sum(y*eta - log(1+exp(eta)))
    if(loglik[t] - loglik[t-1] < tol) return(list(beta=beta,Convergence=cbind(Iteration=(1:t)-1,Loglikelihood=loglik[1:t])))
  }
  stop("The algorithm has not reached convergence")
}

#' MM algorithm for logistic regression
#'
#'
#' @param x  A \verb{n x p} dimensional design matrix. The intercept must be manually included.
#' @param y  A \verb{n} dimensional binary vector.
#' @param tol Threshold at which the algorithm stops.
#' @param beta_start Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients.
#' @param maxiter Maximum number of iteration. If reached, the algorithm raises an error.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
#' 
#' @export
#' 
logit_MM <- function(x, y,  tol = 1e-16, beta_start = NULL, maxiter=10000){

  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")

  loglik <- numeric(maxiter)

  # Initialization (If null, implicitely initialized at beta=0)
  B <- qr(crossprod(x/4,x)) # Bohning matrix (QR version)
  if(is.null(beta_start)) {beta <- solve(B,crossprod(x,y-0.5))}else{beta <- beta_start}

  # Initialization
  eta        <- c(x%*%beta)
  prob       <- 1/(1+exp(-eta))
  w          <- prob*(1-prob)

  # First value of the likelihood
  loglik[1] <- sum(y*eta - log(1+exp(eta)))

  # Iterative procedure
  for(t in 2:maxiter){
    beta       <- beta + solve(B,crossprod(x,y-prob))
    eta        <- c(x%*%beta)
    prob       <- 1/(1+exp(-eta))
    w          <- prob*(1-prob)
    loglik[t]  <- sum(y*eta - log(1+exp(eta)))
    if(loglik[t] - loglik[t-1] < tol) return(list(beta=beta,Convergence=cbind(Iteration=(1:t)-1,Loglikelihood=loglik[1:t])))
  }
  stop("The algorithm has not reached convergence")
}

#' Elastic net algorithm for logistic regression
#'
#' @export
logit_EN <- function(x, y, alpha, lambda, tol = 1e-16, tol_inner = 1e-7, beta_start = NULL, maxiter=10000, maxiter_inner = 10^4, method="ECM"){
  
  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  
  x <- scale(x,scale=FALSE)                                      # Centering
  x <- cbind(1,t(t(x)/apply(x,2,function(x) sqrt(mean(x^2)))))   # Scaling the predictors

  if(method=="ECM") {
  out <- pcd_EM(X,y,alpha = alpha, lambda = lambda, beta_init = beta_start, maxiter=maxiter, tol=tol)
  } else if(method=="MM"){
  out <- pcd_MM(X,y,alpha = alpha, lambda = lambda, beta_init = beta_start, maxiter=maxiter, 
                maxiter_inner=maxiter_inner, tol=tol, tol_inner=tol_inner)
  } else if(method=="NR"){
  out <- pcd_NR(X,y,alpha = alpha, lambda = lambda, beta_init = beta_start, maxiter=maxiter, 
                maxiter_inner=maxiter_inner, tol=tol, tol_inner=tol_inner)
  } else{
    stop("Please provide a valid estimation method")
  }
  list(beta=out$beta, Convergence=cbind(Iteration=(0:(length(out$logpost)-1)),Loglikelihood=out$logpost))
}

#'@export
logpost_R <- function(eta,y,beta,alpha,lambda){
  logpost(eta,y,beta,alpha,lambda)
}

#' Mean-field variational Bayes algorithm for logistic regression
#'
#'
#' @param x  A \verb{n x p} dimensional design matrix. The intercept must be manually included.
#' @param y  A \verb{n} dimensional binary vector.
#' @param tol Threshold at which the algorithm stops.
#' @param prior A list containing the prior mean and prior covariance
#' @param maxiter Maximum number of iteration. If reached, the algorithm raises an error.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
#' 
#' @export
#' 
#'@export
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
    
    if(abs(lowerbound[t] - lowerbound[t-1]) < tol) return(list(mu=mu_vb,Sigma=Sigma_vb, Convergence=cbind(Iteration=(1:t)-1,Lowerbound=lowerbound[1:t]),xi=xi))
  }
  stop("The algorithm has not reached convergence")
}

#' Stochastic variational Bayes algorithm for logistic regression
#'
#'
#' @param x  A \verb{n x p} dimensional design matrix. The intercept must be manually included.
#' @param y  A \verb{n} dimensional binary vector.
#' @param tol Threshold at which the algorithm stops.
#' @param prior A list containing the prior mean and prior covariance
#' @param maxiter Maximum number of iteration. If reached, the algorithm raises an error.
#' @param tau Delay parameter: a number greater than zero
#' @param kappa The forgetting rate: a number between 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
#' 
#' @export
#' 
#'@export
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

