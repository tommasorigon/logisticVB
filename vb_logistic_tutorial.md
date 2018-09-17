

This tutorial reproduces the illustrative simulation study in Section 3 of the paper [Durante and Rigon (2018). *Conditionally conjugate variational Bayes in logistic models*](https://arxiv.org/abs/1711.06999). The core functions of our implementations are made available in the R file [`logistic.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), which can be downloaded from this repository. All the analyses are performed with a **MacBook Air (OS X Sierra, version 10.13.6)**, using an `R` version **3.5.0**. 

As a first step, let us load in memory the file [`logistic.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R) and the `ggplot2` library.

```r
rm(list=ls())
source("logistic.R") # R file that can be downloaded from this github repository
library(ggplot2)     # Plots
```

## Comparison between CAVI and SVI algorithms

We consider an illustrative dataset having a binary responses `y` simulated from a logistic regression with an intercept term and a single continuous covariate `x`. In simulating the response data, we set both the **intercept** and the **slope**, comprising the coefficients vector `beta`, equal to 1. The values of the covariate `x` is instead simulated uniformly in the interval `(-2,2)`.

We approximate the posterior distribution via **Coordinate Ascent Variational Inference** (CAVI) (see `logit_CAVI` function) and **Stochastic Variational Inference** (SVI) (see `logit_SVI` function), considering a moderately diffuse Gaussian prior for the regression coefficients, as outlined in Section 3 of the paper. According to Section 3.2 in the articole, the implementation of SVI requires also the choice of the number of iterations `iter`, the delay `tau` and the forgetting rate `kappa`, which are pre-specified in the code below. 

```r
# True vector of regression coefficients
beta <- c(1, 1)

# Prior hyperparameters
prior <- list(mu = rep(0,2), Sigma = diag(10,2))

# SVI parameter settings
iter    <- 10^4  # Number of iterations
tau     <- 1     # Delay parameter
kappa   <- 0.75  # Forgetting rate parameter
```

Let us first focus on a simple scenario with a small sample size `n` equal to 20.

### Sample size `n = 20`

We conduct an initial simulation with `n = 20`. Let us first simulate the data, consistent with the above discussion.

```r
n <- 20 # Setting the sample size

set.seed(123)      # Set the seed to make this experiment reproducible
x <- runif(n,-2,2) # Generating the covariates
X <- cbind(1,x)    # Design matrix: the intercept is included
y <- rbinom(n,1,prob = plogis(X%*%beta))
```
As anticipated, the **CAVI** algorithm is performed using the `logit_CAVI` function, whereas the **SVI** algorithm is performed using the `logit_SVI` function, both present in the file [`logistic.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R). In both cases, the final solution is obtained rapidly.

```r
set.seed(1010)     # Set the seed to make this experiment reproducible
CAVI_output <- logit_CAVI(X = X, y = y, prior = prior) # CAVI algorithm
SVI_output  <- logit_SVI(X = X, y = y,  prior = prior,  iter = iter, tau = tau, kappa = kappa) # SVI algorithm
```

Finally, we simulate posterior samples from the **Gaussian approximating posteriors arising from CAVI and SVI optimizations**. These samples are useful to compare the approximating distributions from the two algorithms.

```r
set.seed(100)
beta0_CAVI <- rnorm(10^4, CAVI_output$mu[1], sqrt(CAVI_output$Sigma[1,1])) # Posterior distribution of the intercept with CAVI
beta1_CAVI <- rnorm(10^4, CAVI_output$mu[2], sqrt(CAVI_output$Sigma[2,2])) # Posterior distribution of the slope with CAVI
beta0_SVI  <- rnorm(10^4, SVI_output$mu[1], sqrt(SVI_output$Sigma[1,1]))   # Posterior distribution of the intercept with SVI
beta1_SVI  <- rnorm(10^4, SVI_output$mu[2], sqrt(SVI_output$Sigma[2,2]))   # Posterior distribution of the slope with SVI

data_plot <- data.frame(Posterior = c(beta0_CAVI,beta1_CAVI,beta0_SVI,beta1_SVI), beta = rep(rep(c("Intercept","Slope"),each=10^4),2), Algorithm = rep(c("CAVI","SVI"),each=2*10^4), Sample_size = n)
```

### Sample size `n = 100`, `n = 1000` and `n = 10000`.

To study empirically the behavior of the approximating distributions from CAVI and SVI as the sample sizes `n` increases, let us repeat the above simulation for `n = 100`, `n = 1000` and `n = 10000`. The code is exactly the same as the one considered for `n = 20`, with the only difference that the sample size `n` is progressively increased.

```r
nn <- c(100,1000,10000) # Setting the sample size

for(n in nn){
  set.seed(123)      # Set the seed to make this experiment reproducible
  x <- runif(n,-2,2) # Generating the covariate space
  X <- cbind(1,x)    # Design matrix: the intercept is included
  y <- rbinom(n,1,prob = plogis(X%*%beta))

  set.seed(1010)     # Set the seed to make this experiment reproducible
  CAVI_output <- logit_CAVI(X = X, y = y, prior = prior) # CAVI algorithm
  SVI_output  <- logit_SVI(X = X, y = y,  prior = prior,  iter = iter, tau = tau, kappa = kappa) # SVI algorithm

  set.seed(100)
  beta0_CAVI <- rnorm(10^4, CAVI_output$mu[1], sqrt(CAVI_output$Sigma[1,1])) # Posterior distribution of the intercept with CAVI
  beta1_CAVI <- rnorm(10^4, CAVI_output$mu[2], sqrt(CAVI_output$Sigma[2,2])) # Posterior distribution of the slope with CAVI
  beta0_SVI  <- rnorm(10^4, SVI_output$mu[1], sqrt(SVI_output$Sigma[1,1]))   # Posterior distribution of the intercept with SVI
  beta1_SVI  <- rnorm(10^4, SVI_output$mu[2], sqrt(SVI_output$Sigma[2,2]))   # Posterior distribution of the slope with SVI

  data_plot <- rbind(data_plot,data.frame(Posterior = c(beta0_CAVI,beta1_CAVI,beta0_SVI,beta1_SVI), beta = rep(rep(c("Intercept","Slope"),each=10^4),2), Algorithm =    rep(c("CAVI","SVI"),each=2*10^4), Sample_size = n))
}
```

### Results

The following code reproduces **Figure 1** in the paper. 

```r
ggplot(data=data_plot, aes(x = as.factor(Sample_size), y = Posterior, fill=Algorithm)) + facet_grid(~beta) + geom_boxplot(alpha=0.7) + theme_bw() + scale_fill_grey() + geom_hline(yintercept=1, linetype="dotted") + xlab("Sample size") + ylab("Regression Coefficient")
ggsave("img/CAVI_vs_SVI.png", width=9,height=4)
ggsave("img/CAVI_vs_SVI.pdf", width=9,height=4)
```

![](https://raw.githubusercontent.com/tommasorigon/logisticVB/master/img/CAVI_vs_SVI.png)

As is clear from the above figure, although SVI relies on noisy gradients, the final approximations are similar to the optimal solutions produced by CAVI. Moreover, these approximate posteriors increasingly shrink around the true coefficients as `n`, thus supporting Corollary 1 in the paper, while also highlighting how the entire approximated posterior, and not only its expectation, concentrates around the truth. However, especially for large sample sizes, there might be some differences between the CAVI and the SVI algorithm.
