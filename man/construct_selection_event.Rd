\name{construct_selection_event}
\alias{construct_selection_event}
\title{
  Construct selection event
}
\description{
  This function contructs the selection event by computing c_k, d_k and e_k which are the constants in the quadratic inequalities which characterize the model selection event. The function is used internally by the function solve_selection_event, which returns the intervals of the OLS estimator where the selection event takes place. 
}
\usage{
construct_selection_event(a,b,R_M_k,kappa_M_k,R_M_phat,kappa_M_phat)
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
}
\value{
  \item{c_k}{
    Constant c_k in the quadratic inequality c_k*Z^2+d_k*Z+e_k>=0 which characterizes the model selection event of the selected model compared to model k (see Lemma 1 for details)
  }
  \item{d_k}{
    Constant d_k in the quadratic inequality c_k*Z^2+d_k*Z+e_k>=0 which characterizes the model selection event of the selected model compared to model k (see Lemma 1 for details)
  }
  \item{e_k}{
    Constant e_k in the quadratic inequality c_k*Z^2+d_k*Z+e_k>=0 which characterizes the model selection event of the selected model compared to model k (see Lemma 1 for details)
  }
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
