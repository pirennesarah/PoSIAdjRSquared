\name{fit_all_subset_linear_models}
\alias{fit_all_subset_linear_models}
\title{Fit all subset linear models}
\description{
  This function fits all possible combinations of linear models and returns the selected model based on adjusted R^2.
}
\usage{
 fit_all_subset_linear_models(y, X, intercept)
}
\arguments{
  \item{y}{
    Response vector of type "matrix"  and dimension nx1
  }
  \item{X}{
    Design matrix of type "matrix"  and dimension nxp
  }
  \item{intercept}{
    Logical value: TRUE if fitted models should contain intercept, FALSE if not
  }
}
\value{
  \item{k}{
    Index set included in model k
  }
  \item{best_model}{
    The selected model fit (lm object)
  }
  \item{phat}{
    Index set included in the selected model
  }
  \item{X_M_phat}{
    The design matrix in the selected model
  }
  \item{best_adj_r_squared}{
    The adjusted R^2 value of the selected model
  }
  \item{R_M_phat}{
    The orthogonal projection matrix of the selected model
  }
  \item{kappa_M_phat}{
    Adjustment factor for model complexity kappa of the selected model
  }
  \item{R_M_k}{
    The orthogonal projection matrix of model k}
  \item{kappa_M_k}{
    Adjustment factor for model complexity kappa of model k
  }
}
\references{
  Pirenne, S. and Claeskens, G. (2024). Exact Post-Selection Inference for Adjusted R Squared.
}
\examples{
  # Generate data
  Data <- datagen.norm(seed = 7, n = 100, p = 10, rho = 0, beta_vec = c(1,0.5,0,0.5,0,0,0,0,0,0))
  X <- Data$X
  y <- Data$y

  # Select model
  fit_all_subset_linear_models(y, X, intercept=FALSE)
}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
