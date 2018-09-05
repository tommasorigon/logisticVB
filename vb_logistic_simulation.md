

This document aim to reproduce the results of the paper [Durante and Rigon (2017)](https://arxiv.org/abs/1711.06999). Our implementations are made available through the R functions [`functions.R`](), which can be downloaded from this repository.

All the analyses are performed with a **MacBook Air (OS X Sierra, version 10.13.1)**, using a R version **3.5.0**. Notice that the matrix decompositions involved in this code might differ across operating systems. 

Before starting the analysis, we load in memory some useful libraries, including the `logistic` library.


```r
rm(list=ls())
library(logistic)# Library that can be downloaded from this github repository
library(knitr)   # To produce "nice" table
library(ggplot2) # Plots
library(reshape2)# Reshaping data, used together with ggplot2
```

## Comparison between CAVI and SVI algorithms

We consider a simulated dataset having a binary output `y` and a single continous covariate `x`. We obtain the posterior distribution using both the Coordinate Ascent Variational Inference (CAVI) algorithm (`logit_VB` function of the `logistic` package) and the Stochastic Variational Inference (SVI) algorithm (`logit_SVI` function of the `logistic` package).

### Number of observations `n = 100`

```r
# True vector of coefficients
beta <- c(1, 0.5)


X <- cbind(1,1) # Design matrix
```

