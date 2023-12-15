\name{datagen.norm.intercept}
\alias{datagen.norm.intercept}
\title{Data generation normal with intercept}
\description{
  Function to generate data according to the linear model of the form Y = X*beta + epsilon where the noise epsilon follows a standard normal distribution and the first column of X consists of 1's such that an intercept is included in the model.
}
\usage{
  datagen.norm.intercept(seed, n, p, rho, beta_vec)
}
\arguments{
  \item{seed}{
    Integer for seed
  }
  \item{n}{
    Integer for sample size
  }
  \item{p}{
    Integer for number of variables in the design matrix
  }
  \item{rho}{
    Integer for correlation between variables in the design matrix
  }
  \item{beta_vec}{
    True regression coefficient vector of length p
  }
}
\value{
 \item{X}{
    Design matrix of type "matrix"  and dimension nxp
 }
 \item{y}{
    Response vector of type "matrix"  and dimension nx1
 }
 \item{true_y}{
    True response vector, i.e. without the noise, of type "matrix" and dimension nx1
 }
}
\examples{
  datagen.norm.intercept(seed = 7, n = 100, p = 10, rho = 0, beta_vec = c(1,0.5,0,0.5,0,0,0,0,0,0))
}
\keyword{datagen}
\keyword{models}
\keyword{regression}
