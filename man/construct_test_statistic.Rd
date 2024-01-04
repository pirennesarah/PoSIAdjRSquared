\name{construct_test_statistic}
\alias{construct_test_statistic}
\title{Construct test statistic}
\description{
  This function constructs the OLS estimator of the j'th selected coefficient in the selected model. The functions also returns some useful vectors for post-selection inference (a and b).
}
\usage{
 construct_test_statistic(j, X_M_phat, y, phat, Sigma, intercept)
}
\arguments{
  \item{j}{
    The index of type "integer" of the regression coefficient
  }
  \item{X_M_phat}{
    The design matrix in the selected model
  }
  \item{y}{
    Response vector of type "matrix"  and dimension nx1
  }
  \item{phat}{
    Index set included in the selected model
  }
  \item{Sigma}{
    The variance covariance matrix of dimension nxn of the error in the model
  }
  \item{intercept}{
    Logical value: TRUE if the selected model contains an intercept, FALSE if not
  }
}
\value{
  \item{etaj}{
    Vector of type "matrix" and dimension nx1: useful in orthogonal decomposition of y (see Lemma 1 for details)
  }
  \item{etajTy}{
    The OLS estimator of the j'th selected coefficient in the selected model of type "matrix" and dimension 1x1
  }
  \item{a}{
    Residual vector of type "matrix" and dimension nx1 (see Lemma 1 for details)
  }
  \item{b}{
    Vector of type "matrix" and dimension nx1: useful in orthogonal decomposition of y (see Lemma 1 for details)
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

  # Estimate Sigma from residuals of full model
  full_model <- lm(y ~ 0 + X)
  sigma_hat <- sd(resid(full_model))
  Sigma <- diag(n)*(sigma_hat)^2

  # Construct test statistic
  construct_test_statistic(j = 5, X_M_phat, y, phat, Sigma, intercept=FALSE)
}
\keyword{htest}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
\concept{post-selection inference}
