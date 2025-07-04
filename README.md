# PoSIAdjRSquared

The PoSIAdjRSquared package allows users to calculate p-values and confidence intervals 
for regression coefficients after they have been selected by adjusted R squared in linear models. 
The p-values and confidence intervals are valid after model selection with the same data. 
This allows the user to use all data for both model selection and inference without losing control over the type I error rate. 
The provided tests are more powerful than data splitting, which bases inference on less data since it discards all information used for selection.

## Example

This is a basic example which shows you how to calculate post-selection
p-values and confidence intervals for some generated data. The code is
similarly applicable to real data.

``` r
library(PoSIAdjRSquared)

# Generate data
Data <- datagen.norm(seed = 7, n = 100, p = 10, rho = 0, beta_vec = c(1,0.5,0,0.5,0,0,0,0,0,0))
	X <- Data$X
  	y <- Data$y

selective_inference(y, X, intercept=FALSE, model_set = "fit_all_subset_linear_models", alpha=0.05, confidence_interval=TRUE)

# $summary
# Variable Coefficient  Std_Error     P_Value    CI_Lower  CI_Upper
# 1        1   1.2550561 0.10171810 0.000000000  1.00651134 1.4543199
# 2        2   0.3710123 0.10468937 0.322730857 -0.04019441 0.5530123
# 3        4   0.3291952 0.09248687 0.001782217  0.12471371 0.5104508
# 4        5  -0.1234033 0.10508632 0.841042743 -0.23945366 0.1111173
# 5        8   0.1358987 0.09710654 0.548071861 -0.08495766 0.3009992
# 6       10   0.1196511 0.09917412 0.997850742 -0.10880178 0.2263758
```

## Reference
Pirenne, S. and Claeskens, G. (2024). Exact post-selection inference for adjusted R squared selection. Statistics & Probability Letters, 211(110133):1-9. https://doi.org/10.1016/j.spl.2024.110133
