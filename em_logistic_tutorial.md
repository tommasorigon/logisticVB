This tutorial provides quantitative support to the theoretical results in the Appendix of the paper [Durante and Rigon (2018). *Conditionally conjugate variational Bayes in logistic models*](https://arxiv.org/abs/1711.06999). Here the focus is on maximum likelihood estimation of the logistic coefficients and not on Bayesian inference. The core functions of our implementations are made available in the R file [`logistic.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), which can be downloaded from this repository. All the analyses are performed with a **MacBook Air (OS X Sierra, version 10.13.6)**, using a R version **3.5.0**. 

As a first step, let us load in memory the file [`logistic.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), as well as the `ggplot2` and the `knitr` libraries.

```r
rm(list=ls())
source("logistic.R") # R file that can be downloaded from this github repository
library(ggplot2)     # Plots
library(knitr)       # To produce tables
```

## Failure of the Newton-Raphson algorithm in logistic regression

The main aim of this illustrative simulation is to show that the **Newton-Raphson** algorithm for logistic regression, sometimes called also *Fisher scoring*, can fail even when the maximum likelihood estimate (MLE) is well defined. To clarify, this problem is not related to the so-called *separability issue*. In fact, in this case the MLE simply does not exists. To quantitatively evaluate the theoretical results on convergence rates discussed in **Proposition 1** of the paper, we additionally implement the **EM** based on the Pòlya-Gamma data augmentation (see `logit_EM` function) and the **MM** of [Böhning and Lindsay (1988)](https://link.springer.com/article/10.1007/BF00049423) (see `logit_MM` function). 

The dataset considered in this analysis has been suggested in [this blog post](http://www.win-vector.com/blog/2012/08/how-robust-is-logistic-regression/) and comprises a binary response `y` along with an intercept term and a single covariate defining the design matrix `X`. Let us create this dataset below. 

```r
y <- c(rep(0,50),1,rep(0,50),0,rep(0,5),rep(1,10))        # Binary outcomes
X <- cbind(1,c(rep(0,50),0,rep(0.001,50),100,rep(-1,15))) # Design matrix
```

To perform maximum likelihood estimation, we initialize the three algorithms in `(0,0)`, which is the typical choice in standard statistical packages. This can be done by setting `beta_start = c(0,0)`. 

```r
fit_NR <- logit_NR(X,y,beta_start=c(0,0))               # Newton-Raphson
fit_EM <- logit_EM(X,y,beta_start=c(0,0))               # EM via Polya-gamma
fit_MM <- logit_MM(X,y,beta_start=c(0,0), maxiter=10^5) # Böhning and Lindsay (1988)
```

For this specific dataset, the maximum value is attained when the intercept and the slope coefficients are `-4.603`and `-5.296`, respectively, whereas the log-likelihood, evaluated in the maximum, is equal to `-15.156`. To study the performance of the three maximization methods let us study the log-likelihood sequence of the first 5 iterations.

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

Since the algorithms are initialized at the same point `(0,0)`, the log-likelihood at iteration `0` is the same for all the routines. Moreover, they also coincide at iteration `1`, because the three methods produce the same updating step. However, after iteration `2` they become different. The log-likelihood sequence of the **Newton-Raphson** algorithm falls at iteration `5`, and the coefficients estimates diverge. Indeed, note that also the default `glm` R command fails to reach convergence. In fact, the algorithm adopted in `glm` is the Newton-Raphson. 

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

As discussed in the Appendix of the paper, both the **MM** and the **EM** algorithms produce a monotone log-likelihood sequence and eventually reach convergence. However, as is clear from the above table, the EM based on Pòlya-Gamma data augmentation seems to reach the maximum more quickly than MM. This is consistent with **Proposition 1** in the article. Indeed, in this specific simulation the MM requires `41396` iterations to reach convergence, whereas the EM only `360`.

## Empirical convergence rate of the MM and EM algorithms

Let us confirm the above findings on the convergence rates of **MM** and **EM** in a more general setting with `n = 10000` statistical units and several covariates simulated uniformly in the interval `(-2,2)`. The binary response data are instead generated from a logistic regression with true coefficients vector `(1, 1, -1, 1,-1, 1, -1)`.

Let us create this simulated dataset.

```r
n <- 10000 # Setting the sample size

beta <- c(1, 1, -1, 1,-1, 1, -1) # True value of the simulated coefficients

set.seed(123) # Set the seed to make this experiment reproducible
X <- cbind(1,runif(n,-2,2),runif(n,-2,2),runif(n,-2,2),runif(n,-2,2),runif(n,-2,2), runif(n,-2,2))    # Design matrix: the intercept is included
y <- rbinom(n,1,prob = plogis(X%*%beta))
```

We now maximize the log-likelihood via EM and MM. As done before, we initialize the coeffients values at `0`.

```r
fit_EM <- logit_EM(X,y) # EM via Polya-gamma
fit_MM <- logit_MM(X,y) # Böhning and Lindsay (1988)
```

For this specific dataset, the log-likelihood evaluated in the maximum is equal to `-3571.163`. The MM algorithm requires `130` iterations to reach convergence whereas the EM algorithm reaches the maximum more rapidly, requiring `59` iterations. To provide a graphical representation of the improved rate of convergence, we display the first `20` values of the log-likelihood as a function of the iterations.

```r
iters <- 20
data_plot <- data.frame(Iteration = c(fit_EM$Convergence[1:iters,1],fit_MM$Convergence[1:iters,1]), Log_likelihood = c(fit_EM$Convergence[1:iters,2],fit_MM$Convergence[1:iters,2]), Algorithm = as.factor(rep(c("EM","MM"), each = iters)))

ggplot(data=data_plot, aes(x = Iteration, y = Log_likelihood, linetype=Algorithm)) + geom_point(size=0.7) + geom_line() + theme_bw() + xlab("Iteration") + ylab("Log-likelihood")
ggsave("img/EM_vs_MM.png", width=9,height=4)
ggsave("img/EM_vs_MM.pdf", width=9,height=4)
```

![](https://raw.githubusercontent.com/tommasorigon/logisticVB/master/img/EM_vs_MM.png)

As is clear from the above figure, the MM algorithm converges more slowly compared to the EM, although both eventually reach the maximum likelihood estimate. This result confirms again **Proposition 1** in the article.
