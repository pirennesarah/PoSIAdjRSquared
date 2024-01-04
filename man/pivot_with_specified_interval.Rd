\name{pivot_with_specified_interval}
\alias{pivot_with_specified_interval}
\title{Pivot with specified interval
}
\description{
  This function returns the value of the cumulative distribution function of a truncated gaussian distribution evaluated in the observed test statistic. Its output is used by the function postselp_value_specified_interval by taking 2*min(output,1-output) for a p-value for a two-sided test. 
}
\usage{
pivot_with_specified_interval(z_interval, etaj, etajTy, tn_mu, tn_sigma)
}
\arguments{
  \item{z_interval}{
    The intervals of type "list" where the OLS estimator gets selected: can be obtained from function "solve_selection_event"
  }
  \item{etaj}{
    Vector of type "matrix" and dimension nx1: useful in orthogonal decomposition of y (see Lemma 1 for details)
  }
  \item{etajTy}{
    The OLS estimator of the j'th selected coefficient in the selected model of type "matrix" and dimension 1x1
  }
  \item{tn_mu}{
    Integer for the mean of the truncated sampling distribution of the test statistic under the null hypothesis: for example, if you want to test beta_j=0, specify 0 for the mean
  }
  \item{tn_sigma}{
    Integer for the variance of the truncated sampling distribution of the test statistic
  }
}
\value{The cumulative distribution function of a truncated gaussian distribution evaluated in the observed test statistic
}
\references{
  Pirenne, S. and Claeskens, G. (2024). Exact Post-Selection Inference for Adjusted R Squared Selection.
}
\examples{
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

  # Post-selection inference for beta_j=0
  tn_sigma <- sqrt((t(etaj)\%*\%Sigma)\%*\%etaj)
  pivot_with_specified_interval(z_interval, etaj, etajTy, tn_mu = 0, tn_sigma)
}
\keyword{htest}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
\concept{post-selection inference}
