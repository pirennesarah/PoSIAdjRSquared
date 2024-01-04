\name{solve_selection_event}
\alias{solve_selection_event}
\title{Solve selection event}
\description{
  This function solves the selection event by calculating the intervals of the OLS estimator where the regression coefficient actually gets selected.
}
\usage{
 solve_selection_event(a,b,R_M_k,kappa_M_k,R_M_phat,kappa_M_phat,k)
}
\arguments{
  \item{a}{
    Residual vector of type "matrix" and dimension nx1 (see Lemma 1 for details)
  }
  \item{b}{
    Vector of type "matrix" and dimension nx1: useful in orthogonal decomposition of y (see Lemma 1 for details)
  }
  \item{R_M_k}{
    The orthogonal projection matrix of model k
  }
  \item{kappa_M_k}{
    Adjustment factor for model complexity kappa of model k
  }
  \item{R_M_phat}{
    The orthogonal projection matrix of the selected model
  }
  \item{kappa_M_phat}{
    Adjustment factor for model complexity kappa of the selected model
  }
  \item{k}{
    Index set included in model k
  }
}
\value{
  \item{intervals_list}{
    The intervals of the OLS estimator for which the inequality in Lemma 1 holds
  }
  \item{z_interval}{
    The intersection of the intervals of the OLS estimator for which the inequality in Lemma 1 holds: post-selection inference is conditioned on those intervals
  }
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

  # Solve selection event
  solve_selection_event(a,b,R_M_k,kappa_M_k,R_M_phat,kappa_M_phat,k)
}
\keyword{htest}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
\concept{post-selection inference}
\note{
  This function will give the error message: "Error in if (D >= 0.001): missing value where TRUE/FALSE needed" if it is run for a coefficient which is not part of the selected model. Only run solve_selection_event for selected indices.
}
