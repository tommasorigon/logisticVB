This document is associated to the paper [Durante and Rigon (2017). *Conditionally conjugate variational Bayes in logistic models*](https://arxiv.org/abs/1711.06999). The core functions of our implementations are made available in the file [`functions.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), which can be downloaded from this repository.

All the analyses are performed with a **MacBook Air (OS X Sierra, version 10.13.6)**, using a R version **3.5.0**. 

Before starting the analysis, we load in memory the file [`functions.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R) and the `ggplot2` library.

```r
rm(list=ls())
source("logistic.R") # R file that can be downloaded from this github repository
library(ggplot2)     # Plots
library(knitr)       # To produce tables
```

## Failure of the Newton-Raphson algorithm in logistic regression

The aim of this paragraph is to show that the Newton-Raphson algorithm, sometimes called Fisher scoring, can fail even if the maximum lieklihood estimate is well defined. We consider the following dataset having a binary output `y` and continous covariate. 

```r
y <- c(rep(0,50),1,rep(0,50),0,rep(0,5),rep(1,10))        # Binary outcomes
X <- cbind(1,c(rep(0,50),0,rep(0.001,50),100,rep(-1,15))) # Design matrix
```

The maximum likelihood estimate exists and it is well defined. However, the Newton-Raphson algorithm **diverges** after few iterations.  We compare the Newton-Raphson algorithm (`logit_NR` function) with two other algorithms: the EM based on the Polya-Gamma data augmentation (`logit_EM` function) and the MM of [Böhning and Lindsay (1988)](https://link.springer.com/article/10.1007/BF00049423) (`logit_MM` function). The MM and the EM are **monotone** in the sense that they theoretically guarantee that the log-likelihood increases at each iteration.

We initialize the three algorithms in (0,0), which is the typical choice in standard statistical packages. This can be done by setting `beta_start = c(0,0)`. 

```r
fit_NR <- logit_NR(X,y,beta_start=c(0,0))               # Newton-Raphson
fit_EM <- logit_EM(X,y,beta_start=c(0,0))               # EM via Polya-gamma
fit_MM <- logit_MM(X,y,beta_start=c(0,0), maxiter=10^5) # Böhning and Lindsay (1988)
```

For this specific dataset, the maximum values is attained when the beta coefficients are equal to `(-4.603,-5.296)` whereas the log-likelihood, evaluated in the maximum, is equal to `-15.156`.

The MM algorithm requires the huge number of  `41'397` iterations to reach convergence. As expected, the EM algorithm reaches the maximum more rapidly, requiring only `361` iterations. This is not surprising, since the EM algorithm has a faster rate of convergence compared to the MM, as illustrated in the paper. This aspect is also illustrated more in depth in the next paragraph.

Conversely, the Newton-Raphson algorithm fails to reach the global maximum. To illustrate the behaviour of these three algoritms, we show in the following table the value of the log-likelihood of the first 5 iterations.


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

Since the algorithms are initialized at the same point `(0,0)`, the log-likelihood at iteration `0` is the same for all the methods. Moreover, they also coincide at iteration `1`, because the three algorithms produce the same updating step. However, from iteration `2` they become different. The log-likelihood sequence of Newton-Raphson algorithm falls at iteration `5`, and the beta coefficients diverge. Conversely, both the MM and the EM algorithms produce a monotonic sequence and eventually reach convergence.

Notice that also the default `glm` R command fails to reach convergence, although it raises a warning. Indeed, the algorithm adopted in `glm` is the Newton-Raphson. 


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

## Empirical convergence rate of the MM and EM algorithms

In this Section we empirically check that the convergence rate of the EM based on the Polya-gamma distribution has a faster convergence rate compared to the MM approach of Böhning. To this extent, we consider the following simulated dataset, with a a sample size `n = 10000` and several covariates.

```r
n <- 10000 # Setting the sample size

beta <- c(1, 1, -1, 1,-1, 1, -1) # True value of the simulated coefficients

set.seed(123) # Set the seed to make this experiment reproducible
X <- cbind(1,runif(n,-2,2),runif(n,-2,2),runif(n,-2,2),runif(n,-2,2),runif(n,-2,2), runif(n,-2,2))    # Design matrix: the intercept is included
y <- rbinom(n,1,prob = plogis(X%*%beta))
```

Then, we fit a logistic regression model using the EM and the MM algorithm. We initialize both the algorithms at `0`, as before.

```r
fit_EM <- logit_EM(X,y) # EM via Polya-gamma
fit_MM <- logit_MM(X,y) # Böhning and Lindsay (1988)
```

For this specific dataset, the log-likelihood evaluated in the maximum is equal to `-3571.163`. The MM algorithm requires `132` iterations to reach convergence whereas the EM algorithm reaches the maximum more rapidly, requiring `60` iterations. To provide a graphical representation of the increased rate of convergence, we display the value of the log-likelihood as a function of the iterations.

```r
iters <- 20
data_plot <- data.frame(Iteration = c(fit_EM$Convergence[1:iters,1],fit_MM$Convergence[1:iters,1]), Log_likelihood = c(fit_EM$Convergence[1:iters,2],fit_MM$Convergence[1:iters,2]), Algorithm = as.factor(rep(c("EM","MM"), each = iters)))

ggplot(data=data_plot, aes(x = Iteration, y = Log_likelihood, linetype=Algorithm)) + geom_point(size=0.7) + geom_line() + theme_bw() + xlab("Iteration") + ylab("Log-likelihood")
ggsave("img/EM_vs_MM.png", width=9,height=4)
ggsave("img/EM_vs_MM.pdf", width=9,height=4)
```

![](https://raw.githubusercontent.com/tommasorigon/logisticVB/master/img/EM_vs_MM.png)

As apparent from the plot above, the MM algorithm reach converges more slowly compared to the EM, although both eventually reach the maximum likelihood estimate.