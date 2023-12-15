# PoSIAdjRSquared

The PoSIAdjRSquared package allows users to calculate p-values and confidence intervals 
for regression coefficients after they have been selected by adjusted R squared in linear models. 
The p-values and confidence intervals are valid after model selection with the same data. 
This allows the user to use all data for both model selection and inference without losing control over the type I error rate. 
The provided tests are more powerful than data splitting, which bases inference on less data since it discards all information used for selection.

## Installation

You can install the development version of PoSIAdjRSquared from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
library(devtools)
install_github("pirennesarah/PoSIAdjRSquared")
```

## Example

This is a basic example which shows you how to calculate post-selection
p-values and confidence intervals for some generated data. The code is
similarly applicable to real data.

``` r
library(PoSIAdjRSquared)

  # Generate data
  n <- 100
  Data <- datagen.norm(seed = 7, n, p = 10, rho = 0, beta_vec = c(1,0.5,0,0.5,0,0,0,0,0,0))
  X <- Data$X
  y <- Data$y

  # Select model
  result <- fit_all_subset_linear_models(y, X, intercept=FALSE)
  phat <- result$phat
  X_M_phat <- result$X_M_phat
  k <- result$k
  R_M_phat <- result$R_M_phat
  kappa_M_phat <- result$kappa_M_phat
  R_M_k <- result$R_M_k
  kappa_M_k <- result$kappa_M_k

  # Estimate Sigma from residuals of full model
  full_model <- lm(y ~ 0 + X)
  sigma_hat <- sd(resid(full_model))
  Sigma <- diag(n)*(sigma_hat)^2

  # Construct test statistic
  Construct_test <- construct_test_statistic(j = 5, X_M_phat, y, phat, Sigma, intercept=FALSE)
  a <- Construct_test$a
  b <- Construct_test$b
  etaj <- Construct_test$etaj
  etajTy <- Construct_test$etajTy

  # Solve selection event
  Solve <- solve_selection_event(a,b,R_M_k,kappa_M_k,R_M_phat,kappa_M_phat,k)
  z_interval <- Solve$z_interval

  # Post-selection p-value for beta_j=0
  tn_sigma <- sqrt((t(etaj)%*%Sigma)%*%etaj)
  postselp_value_specified_interval(z_interval, etaj, etajTy, tn_mu = 0, tn_sigma)
#> [1] 0.8410427
  
  # Post-selection (1-alpha)% confidence interval
  compute_ci_with_specified_interval(z_interval, etaj, etajTy, Sigma, tn_mu = 0, alpha = 0.05)
#> [1] -0.2394537  0.1111173
```

## Reference
Pirenne, S. and Claeskens, G. (2024). Exact Post-Selection Inference for Adjusted R Squared.
