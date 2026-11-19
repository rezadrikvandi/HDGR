# HDGR

HDGR is an R package for high-dimensional group regression with tailored penalties. It implements the Case 3 method described in the following paper.

Drikvandi, R. (2026). Group Regression With Tailored Penalties for High-Dimensional Data With or Without Predefined Groups. *Scandinavian Journal of Statistics*, **53**(3), 1081–1100. [https://doi.org/10.1111/sjos.70072](https://doi.org/10.1111/sjos.70072)

The main function, `HDGR()`, estimates all regression coefficients using a ridge penalty for the group of interest and a lasso penalty for the pooled nuisance predictors. It returns debiased estimates, confidence intervals and two-sided tests for the group of interest.

## Installation

Run the following commands in R. If you need to install `devtools` first, run `install.packages("devtools")`.

```r
devtools::install_github("rezadrikvandi/HDGR")
library(HDGR)
```

## Example 1. Simulated data with known groups

The package includes a simulated dataset with 100 observations and 500 predictors. Loading `testdata` makes the predictor matrix `X`, response `Y`, five group index vectors and true coefficients `true_beta` available by name. The first ten predictors form the group of interest, and the error standard deviation used in the simulation is 1.

Supplying only `G1` gives two groups for fitting: the specified group of interest and one nuisance group containing all remaining predictors.

```r
data("testdata", package = "HDGR")
fit <- HDGR(X, Y, G1 = Index_XG1, sigma = 1)

fit$CI
fit$testpower
```

The table shows the Bonferroni-adjusted 95% confidence intervals and test decisions from one run. A decision of 1 indicates rejection of the hypothesis that the coefficient is zero. The defaults select penalties by cross-validation, so the results can vary with the fold assignment. To reproduce your own run, call `set.seed()` with a chosen integer before fitting.

| Predictor | Lower | Upper | Test decision |
|:----------|------:|------:|--------------:|
| X1  |  1.0974 | 1.8389 | 1 |
| X2  |  0.3540 | 1.1372 | 1 |
| X3  |  0.0122 | 0.7948 | 1 |
| X4  |  0.4568 | 1.2393 | 1 |
| X5  | -0.2350 | 0.6221 | 0 |
| X6  |  0.2707 | 1.0134 | 1 |
| X7  | -0.1424 | 0.6641 | 0 |
| X8  |  0.2188 | 1.0316 | 1 |
| X9  | -0.0914 | 0.6817 | 0 |
| X10 |  0.4020 | 1.0827 | 1 |

To supply all five known groups, pass one label per column of `X` and identify the target label. The following call uses the penalties selected above so that the two fits can be compared directly.

```r
groups <- integer(ncol(X))
groups[Index_XG1] <- 1L
groups[Index_XG2] <- 2L
groups[Index_XG3] <- 3L
groups[Index_XG4] <- 4L
groups[Index_XG5] <- 5L

fit.groups <- HDGR(X, Y, groups = groups, target = 1, sigma = 1,
                   lambda1 = fit$lambda1, lambda2 = fit$lambda2)
```

## Example 2. Riboflavin data with unknown groups

The riboflavin dataset contains gene expression measurements for 4,088 genes in 71 observations. The response is the log riboflavin production rate. Omit both `G1` and `groups` to apply group correlation learning (GCL). Here GCL forms 200 groups and takes the 21 genes with the highest absolute marginal correlation scores as the group of interest.

```r
data("riboflavin", package = "HDGR")
set.seed(1)
fit.gcl <- HDGR(riboflavin$x, riboflavin$y, n_groups = 200)

fit.gcl$G1            # Column indices of the selected genes
fit.gcl$betaG1_hat    # Penalised estimates for these genes
fit.gcl$inference.G1  # Debiased estimates, standard errors, intervals and tests
```

Here `sigma` is estimated from the data. The default intervals and tests for GCL use the fitted group and are unadjusted for multiplicity. Set `adjustment = "bonferroni"` to adjust across the target coefficients. See `?HDGR` for all arguments and return values.
