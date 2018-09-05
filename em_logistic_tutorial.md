

This document aim to reproduce the simulation study of the paper [Durante and Rigon (2017). *Conditionally conjugate variational Bayes in logistic models*](https://arxiv.org/abs/1711.06999). The core functions of our implementations are made available in the file [`functions.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), which can be downloaded from this repository.

All the analyses are performed with a **MacBook Air (OS X Sierra, version 10.13.6)**, using a R version **3.5.0**. 

Before starting the analysis, we load in memory the file [`functions.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R) and the `ggplot2` library.

```r
rm(list=ls())
source("logistic.R") # R file that can be downloaded from this github repository
library(ggplot2)     # Plots
```

## Failure of the Newton-Raphson algorithm in logistic regression for finding MLE

We consider the following simulated dataset having binary output and continous covariate. The dataset is generated as follows

```r
y <- c(rep(0,50),1,rep(0,50),0,rep(0,5),rep(1,10))        # Binary outcomes
X <- cbind(1,c(rep(0,50),0,rep(0.001,50),100,rep(-1,15))) # Design matrix
```

The maximum likelihood estimate exists and it is well defined. However, the Newton-Raphson algorithm diverges after few iterations.  We compare the Newton-Raphson algorithm with two monotone algorithms: the EM based on the Polya-Gamma data augmentation and the MM of [Böhning and Lindsay (1988)](https://link.springer.com/article/10.1007/BF00049423), through the functions `logit_NR`, `logit_EM` and `logit_MM` respectively. We initialize the three algorithms in (0,0), by setting `beta_start = c(0,0)`. 

```r
fit_NR <- logit_NR(X,y,beta_start=c(0,0))               # Newton-Raphson
fit_EM <- logit_EM(X,y,beta_start=c(0,0))               # EM via Polya-gamma
fit_MM <- logit_MM(X,y,beta_start=c(0,0), maxiter=10^5) # Böhning and Lindsay (1988)
```

For this specific dataset, the maximum values is attained when the beta coefficients are equal to `(-4.603,-5.296)` whereas the log-likelihood evaluated in the maximum is equal to `-15.156`.

The MM algorithm requires the incredible number of  `41'397` iterations to reach convergence and therefore  the default maximum number of iterations was increased. As expected, the EM algorithm reaches the maximum more rapidly, requiring only `361` iterations. Conversely, the Newton-Raphson algorithm fails to reach the global maximum. To illustrate the behaviour of these three algoritms, we show in the following table the value of the log-likelihood of the first 5 iterations.


```r
tab <- rbind(Newton_Raphson=fit_NR$Convergence[1:6,2],
             EM=fit_EM$Convergence[1:6,2],
             MM=fit_MM$Convergence[1:6,2])
colnames(tab) <- 0:5
kable(tab)
```

|               |         0|         1|         2|         3|         4|          5|
|:--------------|---------:|---------:|---------:|---------:|---------:|----------:|
|Newton_Raphson | -81.09822| -38.81425| -36.27081| -35.43309| -26.31365| -733.67090|
|EM             | -81.09822| -38.81425| -36.77784| -36.33180| -36.16827|  -36.06429|
|MM             | -81.09822| -38.81425| -37.02854| -36.52533| -36.33067|  -36.23546|

Since the algorithms are initialized at the same point (0,0), the log-likelihood at iteration `0` is the same for all the methods. Moreover, they also coincide at iteration `1`, because the three algorithms produce the same updating step whenever beta is equal to zero. Starting form iteration `2` they become different. The log-likelihood sequence of Newton-Raphson algorithm falls at iteration `5`, and the beta coefficients diverge. Conversely, both the MM and the EM algorithms produce a monotonic sequence and eventually reach convergence.

Notice that also the `glm` R command fails to reach convergence, although it raises a warning. Indeed, the Fisher-scoring algorithm adopted for glm coincide with the Newton-Raphson in the case of the logistic regression. 


```r
coef(glm(y~X[,-1],family="binomial"))
```

```
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
##   (Intercept)       X[, -1] 
## -3.372166e+15 -2.085057e+13
```

## Convergence rate of MM and EM