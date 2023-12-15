\name{compute_ci_with_specified_interval}
\alias{compute_ci_with_specified_interval}
\title{Compute post-selection confidence interval with specified interval}
\description{
  This function inverts a post-selection p-value to a confidence interval.
}
\usage{
 compute_ci_with_specified_interval(z_interval,etaj,etajTy,Sigma,tn_mu,alpha)
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
  \item{Sigma}{
    The variance covariance matrix of dimension nxn of the error in the model
  }
  \item{tn_mu}{
    Integer for the mean of the truncated sampling distribution of the test statistic under the null hypothesis: for example, if you want to test beta_j=0, specify 0 for the mean
  }
  \item{alpha}{
    Integer for the desired significance level of the confidence interval
  }
}
\value{
  \item{ci}{
    The two-sided (1-alpha)\% confidence interval valid after model selection
  }
}
\references{
  Pirenne, S. and Claeskens, G. (2024). Exact Post-Selection Inference for Adjusted R Squared.
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

  # Post-selection confidence interval
  compute_ci_with_specified_interval(z_interval, etaj, etajTy, Sigma, tn_mu = 0, alpha = 0.05)
}
\keyword{htest}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
\concept{post-selection inference}
