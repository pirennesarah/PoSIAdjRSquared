\name{find_root}
\alias{find_root}
\title{Find root}
\description{
  This function is used internally by the function compute_ci_with_specified_interval for inverting post-selection p-values to confidence intervals. 
}
\usage{
find_root(z_interval, etajTy, tn_sigma, y, lb, ub, tol=1e-6)
}
\arguments{
  \item{z_interval}{
    The intervals of type "list" where the OLS estimator gets selected: can be obtained from function "solve_selection_event"
}
 \item{etajTy}{
    The OLS estimator of the j'th selected coefficient in the selected model of type "matrix" and dimension 1x1
} 
  \item{tn_sigma}{
    Integer for the variance of the truncated sampling distribution of the test statistic
} 
  \item{y}{
    For example 1.0-0.5*alpha for finding the lower bound of a (1-alpha)\% confidence interval, and 0.5*alpha for finding the upper bound of a (1-alpha)\% confidence interval
} 
  \item{lb}{
    Lower bound in current iteration
}
  \item{ub}{
    Upper bound in current iteration
}
  \item{tol}{
    Tolerance parameter: default set to 1e-6
  }
}
\value{ 
    Returns confidence interval bound
}
\references{
  Pirenne, S. and Claeskens, G. (2024). Exact Post-Selection Inference for Adjusted R Squared Selection.
}
\keyword{htest}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
\concept{post-selection inference}
