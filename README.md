
# Conditionally conjugate variational Bayes in logistic models


This repository is associated with the paper [Durante and Rigon (2017). *Conditionally conjugate variational Bayes in logistic models*](https://arxiv.org/abs/1711.06999), and aims at providing detailed materials to fully reproduce the simulation study in the paper. 

The core functions of our implementations are made available in the file [`functions.R`](https://github.com/tommasorigon/logisticVB/blob/master/logistic.R), which can be downloaded from this repository. In addition, we provide two tutorial aimed at fully reproducing our results:

- [`vb_logistic_tutorial.md`](https://github.com/tommasorigon/logisticVB/blob/master/vb_logistic_tutorial.md). In this tutorial, we empirically compare the CAVI and the SVI algorithms on a simulated dataset. We show that the variational approximation arising from both the algorithms concentrates around the "true value" as the sample size increases.
- [`em_logistic_tutorial.md`](https://github.com/tommasorigon/logisticVB/blob/master/em_logistic_tutorial.md). This tutorial has two purpouses. In first place, we show that the Newton-Raphson algorithm for finding the maximum likelihood  estimate (MLE) can fail (i.e. it diverges), even though the MLE is well defined. In second place, we empirically compare the speed of convergence of the EM algorithm based on a Polya-gamma data augmentation, with the MM algorithm of [Böhning and Lindsay (1988)](https://link.springer.com/article/10.1007/BF00049423). We show that the EM approach is faster than the MM in terms of number of iterations needed to reach convergence, which is in line the the theoretical findings of our paper.