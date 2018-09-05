

This document aim to reproduce the simulation study of the paper [Durante and Rigon (2017). *Conditionally conjugate variational Bayes in logistic models*](https://arxiv.org/abs/1711.06999). The core functions of our implementations are made available in the file [`functions.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), which can be downloaded from this repository.

All the analyses are performed with a **MacBook Air (OS X Sierra, version 10.13.1)**, using a R version **3.5.0**. 

Before starting the analysis, we load in memory some useful libraries, including the `logistic` library.


```r
rm(list=ls())
source("logistic.R") # R file that can be downloaded from this github repository
library(knitr)       # To produce "nice" table
library(ggplot2)     # Plots
library(reshape2)    # Reshaping data, used together with ggplot2
```

## Comparison between CAVI and SVI algorithms

We consider a simulated dataset having a binary output `y` and a single continous covariate `x`. We obtain the posterior distribution using both the Coordinate Ascent Variational Inference (CAVI) algorithm (`logit_VB` function) and the Stochastic Variational Inference (SVI) algorithm (`logit_SVI` function).

In the following code, we set the unknown regression coefficients `beta`, as well as the prior hyperparameters. Moreover, the SVI algorithm requires the choice of some tuning parameters: the number of iterations (`maxiter`), the delay (`tau`) and the forgetting rate (`kappa`). All these parameter settings are fixed through the simulations. 

```r
# True vector of regression coefficients
beta <- c(1, 0.5)

# Prior hyperparameters
prior <- list(mu = rep(0,2), Sigma = diag(10,2))

# SVI parameter settings
iter    <- 10^4 # Number of iterations
tau     <- 0    # Delay parameter
kappa   <- 0.7  # Forgetting rate parameter
```

We replicate the simulation study for different sample sizes `n`, to assess the performance of the variational approximations to concentrate around the true value `beta` as `n` increases. 

### Number of observations `n = 100`

```r
# True vector of coefficients
beta <- c(1, 0.5)


X <- cbind(1,1) # Design matrix
```

