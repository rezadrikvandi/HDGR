This R package, HDGR, is for the following published paper: 

Drikvandi, R. (2026). Group Regression With Tailored Penalties for High-Dimensional Data With or Without Predefined Groups. *Scandinavian Journal of Statistics*, 53(3), 1081–1100. [https://doi.org/10.1111/sjos.70072](https://doi.org/10.1111/sjos.70072)

You can install this R package using the following two commands in R:

library(devtools)

install_github("rezadrikvandi/HDGR")

and then load it using library(HDGR)

Details about this package can be found in the DESCRIPTION file. The main function, `HDGR()`, estimates all regression coefficients using a ridge penalty for the group of interest and a lasso penalty for the pooled nuisance predictors (main method, Case 3 in the paper). It returns the corresponding debiased estimates, confidence intervals and two-sided tests for the group of interest.

## Example 1. Simulated data with known groups

The package includes a simulated dataset with 100 observations and 500 predictors. Loading `testdata` makes the predictor matrix `X`, response `Y`, five group index vectors and true coefficients `true_beta` available by name. The first ten predictors form the group of interest, and the error standard deviation used in the simulation is 1.

Supplying only `G1` gives two groups for fitting: the specified group of interest and one nuisance group containing all remaining predictors.

```r
data("testdata", package = "HDGR")
fit <- HDGR(X, Y, G1 = Index_XG1, sigma = 1)

fit$CI
fit$testpower
```

Also, to supply all five known groups in this simulated data example, pass one label per column of `X` and identify the target label. The following call uses the penalties selected above so that the two fits can be compared directly.

```r
groups <- integer(ncol(X))
groups[Index_XG1] <- 1L
groups[Index_XG2] <- 2L
groups[Index_XG3] <- 3L
groups[Index_XG4] <- 4L
groups[Index_XG5] <- 5L

fit.groups <- HDGR(X, Y, groups = groups, target = 1, sigma = 1,
                   lambda1 = fit$lambda1, lambda2 = fit$lambda2)
fit.groups
```

## Example 2. Real data with unknown groups

The riboflavin dataset described in the paper contains gene expression measurements for 4,088 genes in 71 observations. The response is the log riboflavin production rate. Omit both `G1` and `groups` to apply group correlation learning (GCL). Here GCL forms 200 groups and takes the 21 genes with the highest absolute marginal correlation scores as the group of interest.

```r
data("riboflavin", package = "HDGR")
set.seed(1)
fit.gcl <- HDGR(riboflavin$x, riboflavin$y, n_groups = 200)

fit.gcl$G1            # Column indices of the selected genes
fit.gcl$betaG1_hat    # Penalised estimates for these genes
fit.gcl$inference.G1  # Debiased estimates, standard errors, intervals and tests
```

Here `sigma` is estimated from the data as in the paper. The default intervals and tests for GCL use the fitted group and are unadjusted for multiplicity. Set `adjustment = "bonferroni"` to adjust across the target coefficients.
