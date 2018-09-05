

This document aim to reproduce the simulation study of the paper [Durante and Rigon (2017). *Conditionally conjugate variational Bayes in logistic models*](https://arxiv.org/abs/1711.06999). The core functions of our implementations are made available in the file [`functions.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), which can be downloaded from this repository.

All the analyses are performed with a **MacBook Air (OS X Sierra, version 10.13.1)**, using a R version **3.5.0**. 

Before starting the analysis, we load in memory some useful libraries, including the `logistic` library.

```r
rm(list=ls())
source("logistic.R") # R file that can be downloaded from this github repository
library(knitr)       # To produce "nice" table
library(ggplot2)     # Plots
```

## Comparison between CAVI and SVI algorithms

We consider a simulated dataset having a binary output `y` and a single continous covariate `x`. We obtain the posterior distribution using both the Coordinate Ascent Variational Inference (CAVI) algorithm (`logit_CAVI` function) and the Stochastic Variational Inference (SVI) algorithm (`logit_SVI` function).

In the following code, we set the unknown regression coefficients `beta`, as well as the prior hyperparameters. Moreover, the SVI algorithm requires the choice of some tuning parameters: the number of iterations (`iter`), the delay (`tau`) and the forgetting rate (`kappa`). All these parameter settings are fixed through the simulations. 

```r
# True vector of regression coefficients
beta <- c(1, 1)

# Prior hyperparameters
prior <- list(mu = rep(0,2), Sigma = diag(10,2))

# SVI parameter settings
iter    <- 10^4 # Number of iterations
tau     <- 0    # Delay parameter
kappa   <- 0.7  # Forgetting rate parameter
```

We replicate the simulation study for different sample sizes `n`, to assess the performance of the variational approximations to concentrate around the true value `beta` as `n` increases. 

### Number of observations `n = 100`

We conduct the simulation for `n = 100`. The covariate `x` is randomly generated taking uniform values over the space (-2,2).

```r
n <- 100 # Setting the sample size

set.seed(123)     # Set the seed to make this experiment reproducible
x <- runif(n,-2,2) # Generating the covariate space
X <- cbind(1,x)    # Design matrix: the intercept is included
y <- rbinom(n,1,prob = plogis(X%*%beta))
```
As anticipated, the CAVI algorithm is performed using the `logit_CAVI` function, whereas the SVI algorithm is performed using the `logit_SVI` function. In both cases, the final solution is obtained quite rapidly.

```r
set.seed(1010)     # Set the seed to make this experiment reproducible
CAVI_output <- logit_CAVI(X = X, y = y, prior = prior) # CAVI algorithm
SVI_output  <- logit_SVI(X = X, y = y,  prior = prior,  iter = iter, tau = tau, kappa = kappa) # SVI algorithm

# Show the estimated mean using CAVI and SVI
kable(data.frame(CAVI_mean = CAVI_output$mu, SVI_mean = SVI_output$mu))
```

In the above Table we show the mean of the variational solution for the CAVI algorithm and the SVI algorithm.

| CAVI_mean|  SVI_mean|
|---------:|---------:|
|  1.312584| 1.3456096|
|  0.949651| 0.9089779|

Finally, we simulate posterior draws from the CAVI variational distribution from the SVI variational distribution. These values will be used to construct the final plot.

```r
set.seed(100)
beta0_CAVI <- rnorm(10^4, CAVI_output$mu[1], sqrt(CAVI_output$Sigma[1,1]))
beta1_CAVI <- rnorm(10^4, CAVI_output$mu[2], sqrt(CAVI_output$Sigma[2,2]))

beta0_SVI <- rnorm(10^4, SVI_output$mu[1], sqrt(SVI_output$Sigma[1,1]))
beta1_SVI <- rnorm(10^4, SVI_output$mu[2], sqrt(SVI_output$Sigma[2,2]))

data_plot <- data.frame(Posterior = c(beta0_CAVI,beta1_CAVI,beta0_SVI,beta1_SVI), beta = rep(rep(c("Intercept","Slope"),each=10^4),2), Algorithm = rep(c("CAVI","SVI"),each=2*10^4), Sample_size = n)
```

### Number of observations `n = 1'000`, `n = 10'000` and `n = 100'000`.

```r
n <- 1000 # Setting the sample size

set.seed(123)     # Set the seed to make this experiment reproducible
x <- runif(n,-2,2) # Generating the covariate space
X <- cbind(1,x)    # Design matrix: the intercept is included
y <- rbinom(n,1,prob = plogis(X%*%beta))

set.seed(1010)     # Set the seed to make this experiment reproducible
CAVI_output <- logit_CAVI(X = X, y = y, prior = prior) # CAVI algorithm
SVI_output  <- logit_SVI(X = X, y = y,  prior = prior,  iter = iter, tau = tau, kappa = kappa) # SVI algorithm

set.seed(100)
beta0_CAVI <- rnorm(10^4, CAVI_output$mu[1], sqrt(CAVI_output$Sigma[1,1]))
beta1_CAVI <- rnorm(10^4, CAVI_output$mu[2], sqrt(CAVI_output$Sigma[2,2]))

beta0_SVI <- rnorm(10^4, SVI_output$mu[1], sqrt(SVI_output$Sigma[1,1]))
beta1_SVI <- rnorm(10^4, SVI_output$mu[2], sqrt(SVI_output$Sigma[2,2]))

data_plot <- rbind(data_plot,data.frame(Posterior = c(beta0_CAVI,beta1_CAVI,beta0_SVI,beta1_SVI), beta = rep(rep(c("Intercept","Slope"),each=10^4),2), Algorithm = rep(c("CAVI","SVI"),each=2*10^4), Sample_size = n))
```

```r
n <- 10000 # Setting the sample size

set.seed(123)     # Set the seed to make this experiment reproducible
x <- runif(n,-2,2) # Generating the covariate space
X <- cbind(1,x)    # Design matrix: the intercept is included
y <- rbinom(n,1,prob = plogis(X%*%beta))

set.seed(1010)     # Set the seed to make this experiment reproducible
CAVI_output <- logit_CAVI(X = X, y = y, prior = prior) # CAVI algorithm
SVI_output  <- logit_SVI(X = X, y = y,  prior = prior,  iter = iter, tau = tau, kappa = kappa) # SVI algorithm

set.seed(100)
beta0_CAVI <- rnorm(10^4, CAVI_output$mu[1], sqrt(CAVI_output$Sigma[1,1]))
beta1_CAVI <- rnorm(10^4, CAVI_output$mu[2], sqrt(CAVI_output$Sigma[2,2]))

beta0_SVI <- rnorm(10^4, SVI_output$mu[1], sqrt(SVI_output$Sigma[1,1]))
beta1_SVI <- rnorm(10^4, SVI_output$mu[2], sqrt(SVI_output$Sigma[2,2]))

data_plot <- rbind(data_plot,data.frame(Posterior = c(beta0_CAVI,beta1_CAVI,beta0_SVI,beta1_SVI), beta = rep(rep(c("Intercept","Slope"),each=10^4),2), Algorithm = rep(c("CAVI","SVI"),each=2*10^4), Sample_size = n))
```

```r
n <- 100000 # Setting the sample size

set.seed(123)     # Set the seed to make this experiment reproducible
x <- runif(n,-2,2) # Generating the covariate space
X <- cbind(1,x)    # Design matrix: the intercept is included
y <- rbinom(n,1,prob = plogis(X%*%beta))

set.seed(1010)     # Set the seed to make this experiment reproducible
CAVI_output <- logit_CAVI(X = X, y = y, prior = prior) # CAVI algorithm
SVI_output  <- logit_SVI(X = X, y = y,  prior = prior,  iter = iter, tau = tau, kappa = kappa) # SVI algorithm

set.seed(100)
beta0_CAVI <- rnorm(10^4, CAVI_output$mu[1], sqrt(CAVI_output$Sigma[1,1]))
beta1_CAVI <- rnorm(10^4, CAVI_output$mu[2], sqrt(CAVI_output$Sigma[2,2]))

beta0_SVI <- rnorm(10^4, SVI_output$mu[1], sqrt(SVI_output$Sigma[1,1]))
beta1_SVI <- rnorm(10^4, SVI_output$mu[2], sqrt(SVI_output$Sigma[2,2]))

data_plot <- rbind(data_plot,data.frame(Posterior = c(beta0_CAVI,beta1_CAVI,beta0_SVI,beta1_SVI), beta = rep(rep(c("Intercept","Slope"),each=10^4),2), Algorithm = rep(c("CAVI","SVI"),each=2*10^4), Sample_size = n))
```



### Final plot

```r
data_plot$Sample_size <- as.factor(data_plot$Sample_size)
ggplot(data=data_plot, aes(x = Sample_size, y = Posterior, fill=Algorithm)) + facet_grid(~beta) + geom_boxplot(alpha=0.7) + theme_bw() + scale_fill_grey() + geom_hline(yintercept=1, linetype="dotted") + xlab("Sample size") + ylab("Regression Coefficient")
```


