# QuantumComputingBayesianGST

The following code is written in R and stan. Code that implements this is written in Matlab.

## Purpose 
Frequentist statistics requires Monte Carlo simulations in order to get error bars on estimates. For large and computationally difficult problems, this might require a cluster. However, this can be avoided with Bayesian statistics, which provides error bars on estimates. Rather than use MLE to perform Gate Set Tomogrpahy, this code uses Bayesian Statistics. This code is initially run using Matlab code (which cannot be posted). The Matlab code runs the R code, which implements Rstan. The results of Rstan are then fed into Matlab to complete the analysis. 
