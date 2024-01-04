\name{f}
\alias{f}
\title{f}
\description{
Function used internally by compute_ci_with_specified_interval for calculating valid confidence intervals post-selection. 
}
\arguments{
  \item{z_interval}{
    The intervals of type "list" where the OLS estimator gets selected: can be obtained from function "solve_selection_event"
}
  \item{etajTy}{
    The OLS estimator of the j'th selected coefficient in the selected model of type "matrix" and dimension 1x1
}
  \item{mu}{
    Integer for the mean of the truncated sampling distribution of the test statistic (updated iteratively in compute_ci_with_specified_interval)
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
\keyword{htest}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
\concept{post-selection inference}
