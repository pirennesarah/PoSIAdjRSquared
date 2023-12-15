\name{construct_adj_r_squared}
\alias{construct_adj_r_squared}
\title{Construct adjusted R squared}
\usage{
construct_adj_r_squared(X, k, y, n, intercept = c(TRUE, FALSE), sst)
}
\arguments{
\item{X}{Design matrix of type "matrix"  and dimension nxp}

\item{k}{Index set included in model k}

\item{y}{Response vector of type "matrix"  and dimension nx1}

\item{n}{An integer for the sample size}

\item{intercept}{Logical value: TRUE if fitted models should contain intercept, FALSE if not}

\item{sst}{An integer for the total sum of squares}
}
\value{
  \item{X_M_k}{The design matrix of model k}
  \item{P_M_k}{The projection matrix of model k}
  \item{R_M_k}{The orthogonal projection matrix of model k}
  \item{kappa_M_k}{Adjustment factor for model complexity kappa of model k}
  \item{adj_r_squared}{The adjusted R squared value of model k}
}
\description{
This function computes the adjusted R squared and returns some useful matrices from this computation.
}
\references{
  Pirenne, S. and Claeskens, G. (2024). Exact Post-Selection Inference for Adjusted R Squared.
}
\examples{
# Generate data
n <- 100
k <- 1:10
Data <- datagen.norm(seed = 7, n, p = 10, rho = 0, beta_vec = c(1,0.5,0,0.5,0,0,0,0,0,0))
X <- Data$X
y <- Data$y
sst <- sum((y-mean(y))^2)

construct_adj_r_squared(X, k, y, n, intercept=FALSE, sst)
}
\keyword{htest}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
\concept{post-selection inference}
