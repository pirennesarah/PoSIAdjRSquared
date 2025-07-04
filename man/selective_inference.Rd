\name{selective_inference}
\alias{selective_inference}
\title{Selective inference}
\description{
  This function conducts selective inference for selected models by adjusted R squared. Selective p-values (and optionally, selective confidence intervals) are calculated.
}
\usage{
  selective_inference(y, X, intercept, model_set, alpha, confidence_interval, size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    Response vector of type "matrix" and dimension nx1
  }
  \item{X}{
    The full model design matrix
  }
  \item{intercept}{
    Logical value: TRUE if the selected model contains an intercept, FALSE if not
  }
  \item{model_set}{
    Can take values "fit_all_subset_linear_models" and "fit_specified_size_subset_linear_models". For selection out of all possible subset of the full model, choose option "fit_all_subset_linear_models". For selection from models of a specified size, choose option "fit_specified_size_subset_linear_models".
  }
  \item{alpha}{
    The significance level: number between 0 and 1. For example, for 95\% confidence intervals, choose alpha=0.05.
  }
  \item{confidence_interval}{
    Logical value: TRUE if selective confidence intervals should be calculated, FALSE if not
  }
  \item{size}{
    Integer for specified size of selected model (number of coefficients). Use if model_set is "fit_specified_size_subset_linear_models".
  }
}
\value{
  \item{selected_model}{
    The indices of the selected columns of X.
  }
  \item{coefficients}{
    The estimated coefficients in the selected model.
  }
  \item{standard_errors}{
    The standard errors of the coefficients in the selected model.
  }
  \item{p_values}{
    The selective p-values for tests whether the selected coefficients are 0 with two-sided alternatives.
  }
  \item{confidence_intervals}{
    The selective confidence intervals for the coefficients in the selected model.
  }
  \item{summary}{
    A summary table.
  }
}
\references{
  Pirenne, S. and Claeskens, G. (2024). Exact Post-Selection Inference for Adjusted R Squared.
}
\examples{
  # Generate example data
  Data <- datagen.norm(seed = 7, n = 100, p = 10, rho = 0, beta_vec = c(1,0.5,0,0.5,0,0,0,0,0,0))
  	X <- Data$X
    	y <- Data$y

  selective_inference(y, X, intercept=FALSE, model_set = "fit_all_subset_linear_models",
                      alpha=0.05, confidence_interval=TRUE)

  # $summary
  # Variable Coefficient  Std_Error     P_Value    CI_Lower  CI_Upper
  # 1        1   1.2550561 0.10171810 0.000000000  1.00651134 1.4543199
  # 2        2   0.3710123 0.10468937 0.322730857 -0.04019441 0.5530123
  # 3        4   0.3291952 0.09248687 0.001782217  0.12471371 0.5104508
  # 4        5  -0.1234033 0.10508632 0.841042743 -0.23945366 0.1111173
  # 5        8   0.1358987 0.09710654 0.548071861 -0.08495766 0.3009992
  # 6       10   0.1196511 0.09917412 0.997850742 -0.10880178 0.2263758
}
\keyword{models}
\keyword{regression}
\concept{model selection}
\concept{adjusted R squared}
\concept{post-selection inference}
\author{
 Sarah Pirenne and Gerda Claeskens.
}
