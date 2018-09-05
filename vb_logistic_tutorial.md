

This document aim to reproduce the simulation study of the paper [Durante and Rigon (2017). *Conditionally conjugate variational Bayes in logistic models*](https://arxiv.org/abs/1711.06999). The core functions of our implementations are made available in the file [`functions.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), which can be downloaded from this repository.

All the analyses are performed with a **MacBook Air (OS X Sierra, version 10.13.1)**, using a R version **3.5.0**. 

Before starting the analysis, we load in memory some useful libraries, including the `logistic` library.


```r
rm(list=ls())
source("logistic.R") # R file that can be downloaded from this github repository
library(knitr)      # To produce "nice" table
library(ggplot2)    # Plots
library(reshape2)   # Reshaping data, used together with ggplot2
```

## Comparison between CAVI and SVI algorithms

We consider a simulated dataset having a binary output `y` and a single continous covariate `x`. We obtain the posterior distribution using both the Coordinate Ascent Variational Inference (CAVI) algorithm (`logit_VB` function of the `logistic` package) and the Stochastic Variational Inference (SVI) algorithm (`logit_SVI` function of the `logistic` package).

### Number of observations `n = 100`

```r
# True vector of coefficients
beta <- c(1, 0.5)


X <- cbind(1,1) # Design matrix
```

