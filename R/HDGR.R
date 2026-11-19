#' Group Regression With Tailored Penalties for High-Dimensional Data With or Without Predefined Groups
#'
#' Fits the main method (Case 3) with a ridge penalty for the group of
#' interest and a lasso penalty for the pooled nuisance variables
#'
#' @param X Numeric design matrix or numeric data frame with n observations
#'   and p predictors. The centred target columns must have full column rank.
#' @param Y Numeric response vector of length n. Both X and Y are centred internally.
#' @param G1 Optional column indices or names of the predictors of interest.
#'   Supplying only G1 gives two groups for fitting: G1 and one nuisance group containing all remaining predictors.
#' @param groups Optional vector of p group labels in the column order of X.
#'   Supply either G1 or groups. If both are NULL, group correlation learning
#'   (GCL) determines the groups.
#' @param target Label of the group of interest when groups is supplied.
#'   All other groups are pooled for fitting. Defaults to 1.
#' @param n_groups Number of groups to form by GCL, between 2 and p. Group 1
#'   contains the largest absolute marginal correlation scores. Ties are
#'   resolved in the original column order.
#' @param lambda1 Positive ridge penalty or vector of candidates for
#'   cross-validation. NULL uses the default grid for the grouping method.
#' @param lambda2 Positive lasso penalty or vector of candidates on the glmnet
#'   scale. The transformed nuisance predictors are standardised for fitting.
#' @param nfolds Number of cross-validation folds. Defaults to 10. Estimation
#'   of sigma requires at least three folds.
#' @param foldid Optional vector of n integer fold labels starting at 1.
#'   All penalty pairs use the same folds and the full-data grouping and centring.
#' @param sigma Error standard deviation. If NULL and needed for inference,
#'   a cross-validated lasso fit estimates it as the square root of the residual
#'   sum of squares divided by n minus the number of nonzero slopes.
#' @param alpha Significance level for confidence intervals and two-sided tests.
#' @param inference Use "auto" (default) for Wald intervals and tests based on
#'   the debiased G1 estimates, or "none" to skip them. "working.normal" is
#'   also accepted and gives the same calculation as "auto".
#' @param adjustment Multiplicity adjustment for the G1 intervals and tests.
#'   NULL uses "bonferroni" for supplied groups and "none" for GCL.
#' @param nuisance.inference Logical. If TRUE, also test nuisance coefficients
#'   using the lasso projection method in the hdi package.
#' @param nuisance.sigma Optional error standard deviation for nuisance tests.
#'   Defaults to the supplied or estimated sigma.
#'
#' @return A list including G1 and G2d (target and nuisance column indices),
#'   groups (one label per predictor), betaG1_hat and betaG2d_hat (penalised
#'   estimates), and betaG1_debiased (debiased target estimates).
#'   CI gives the lower and upper confidence limits in G1 order. For each
#'   target coefficient, testpower is 1 if the test rejects a zero coefficient
#'   and 0 otherwise. The inference.G1 table also gives standard errors, test
#'   statistics, and raw and adjusted p-values. These three inference outputs
#'   are NULL when inference is "none". The list also includes the selected
#'   penalties, cross-validation results, intercept, coefficients in original
#'   column order, fitted values, sigma, and any requested nuisance inference.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 60), nrow = 40)
#' Y <- 2 * X[, 1] + X[, 2] + rnorm(40)
#'
#' # Only the group of interest is known.
#' fit <- HDGR(X, Y, G1 = 1:3, sigma = 1, lambda1 = 0.05, lambda2 = 0.1)
#' fit$CI
#' fit$testpower
#'
#' # All group labels are known. Group 1 contains the same three predictors.
#' labels <- rep(1:4, times = c(3, 19, 19, 19))
#' fit.groups <- HDGR(X, Y, groups = labels, target = 1, sigma = 1,
#'                    lambda1 = 0.05, lambda2 = 0.1)
#'
#' # Groups are unknown. GCL forms 20 groups of three predictors each.
#' fit.gcl <- HDGR(X, Y, n_groups = 20, sigma = 1,
#'                 lambda1 = 0.05, lambda2 = 0.1)
#' fit.gcl$G1
#' fit.gcl$inference.G1
#' @export
HDGR <- function(X, Y, G1 = NULL, groups = NULL, target = 1,
                 n_groups = min(200L, ncol(X)), lambda1 = NULL,
                 lambda2 = 10^seq(-0.7, -1.1, by = -0.05),
                 nfolds = 10L, foldid = NULL, sigma = NULL, alpha = 0.05,
                 inference = c("auto", "working.normal", "none"),
                 adjustment = NULL, nuisance.inference = FALSE,
                 nuisance.sigma = NULL) {
  inference <- match.arg(inference)
  positive <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  }
  integer.in <- function(x, lower, upper) {
    is.numeric(x) && length(x) > 0L &&
      all(is.finite(x) & x == floor(x) & x >= lower & x <= upper)
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Please install glmnet.")
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if (!is.numeric(X) || any(!is.finite(X)) || n < 3L || p < 2L)
    stop("X must be finite and numeric, with at least 3 rows and 2 columns.")
  if (!is.numeric(Y) || length(Y) != n || any(!is.finite(Y)))
    stop("Supply one finite numeric response per observation.")
  Y <- as.numeric(Y)
  nm <- colnames(X)
  if (is.null(nm)) nm <- paste0("X", seq_len(p))
  if (anyNA(nm) || anyDuplicated(nm) || any(!nzchar(nm)))
    stop("Column names must be nonempty and unique.")
  colnames(X) <- nm
  xmean <- colMeans(X)
  ymean <- mean(Y)
  X <- sweep(X, 2L, xmean)
  Y <- Y - ymean
  norms <- sqrt(colSums(X^2))
  if (any(norms == 0) || all(Y == 0)) stop("Remove constant predictors or response.")
  if (!positive(alpha) || alpha >= 1) stop("alpha must lie between 0 and 1.")
  if (!is.null(sigma) && !positive(sigma)) stop("sigma must be positive.")
  if (!is.null(nuisance.sigma) && !positive(nuisance.sigma))
    stop("nuisance.sigma must be positive.")
  if (!is.logical(nuisance.inference) || length(nuisance.inference) != 1L ||
      is.na(nuisance.inference)) stop("nuisance.inference must be TRUE or FALSE.")
  if (!is.null(G1) && !is.null(groups)) stop("Supply either G1 or groups.")
  learned <- is.null(G1) && is.null(groups)
  scores <- NULL
  if (learned) {
    if (length(n_groups) != 1L || !integer.in(n_groups, 2, p))
      stop("n_groups must be an integer between 2 and ncol(X).")
    scores <- stats::setNames(abs(drop(crossprod(X, Y))) / norms, nm)
    ranked <- order(scores, seq_len(p))
    sizes <- diff(floor((0:n_groups) * as.double(p) / n_groups))
    groups <- integer(p)
    groups[ranked] <- rep(n_groups:1L, sizes)
    target <- 1L
    G1 <- which(groups == target)
  } else if (!is.null(groups)) {
    if (length(groups) != p || anyNA(groups) || length(target) != 1L || anyNA(target))
      stop("Supply one nonmissing group label per column and one target label.")
    G1 <- which(as.character(groups) == as.character(target))
  } else {
    if (is.character(G1)) G1 <- match(G1, nm)
    if (!integer.in(G1, 1, p) || anyDuplicated(G1))
      stop("G1 must contain distinct column positions or names.")
    groups <- rep(2L, p)
    groups[G1] <- 1L
    target <- 1L
  }
  q <- length(G1)
  if (!q || q >= p) stop("Both the target and nuisance blocks must be nonempty.")
  nuisance <- setdiff(seq_len(p), G1)
  if (q >= n || qr(X[, G1, drop = FALSE])$rank < q)
    stop("The centred G1 design must have full column rank for debiasing.")
  if (is.null(lambda1))
    lambda1 <- 10^seq(-1, if (learned) -1.5 else -1.4, by = -0.05)
  valid.grid <- function(x) is.numeric(x) && length(x) > 0L &&
    all(is.finite(x) & x > 0)
  if (!valid.grid(lambda1) || !valid.grid(lambda2))
    stop("Penalty grids must contain finite positive values.")
  lambda1 <- unique(lambda1)
  lambda2 <- unique(lambda2)
  if (is.null(adjustment)) adjustment <- if (learned) "none" else "bonferroni"
  adjustment <- match.arg(adjustment, c("none", "bonferroni"))
  target.tests <- inference != "none"
  estimate.sigma <- is.null(sigma) &&
    (target.tests || (nuisance.inference && is.null(nuisance.sigma)))
  tune <- length(lambda1) * length(lambda2) > 1L
  if (tune || estimate.sigma) {
    if (is.null(foldid)) {
      if (length(nfolds) != 1L || !integer.in(nfolds, 2, Inf))
        stop("nfolds must be an integer of at least 2.")
      foldid <- sample(rep_len(seq_len(min(nfolds, n)), n))
    }
    if (length(foldid) != n || !integer.in(foldid, 1, n) ||
        !identical(sort(unique(as.integer(foldid))), seq_len(max(foldid))) ||
        max(foldid) < 2 || any(n - tabulate(foldid) < 2))
      stop("Invalid folds: use consecutive labels from 1 with at least two training rows.")
    if (estimate.sigma && max(foldid) < 3) stop("CV variance estimation needs >= 3 folds.")
  }

  # SVD applies (I-P_G1)^(1/2) without forming an n by n matrix.
  prepare <- function(rows) {
    A <- X[rows, G1, drop = FALSE]
    Z <- X[rows, nuisance, drop = FALSE]
    y <- Y[rows]
    s <- svd(A)
    list(Z = Z, y = y, U = s$u, V = s$v, d = s$d,
         Uy = drop(crossprod(s$u, y)), UZ = crossprod(s$u, Z))
  }
  fit.path <- function(a, l1, l2, keep = FALSE) {
    penalty <- length(a$y) * l1
    den <- a$d^2 + penalty
    shrink <- 1 - sqrt(penalty / den)
    Zstar <- a$Z - a$U %*% (shrink * a$UZ)
    ystar <- a$y - drop(a$U %*% (shrink * a$Uy))
    design <- Zstar
    # glmnet requires a matrix with at least two columns.
    if (ncol(design) == 1L) design <- cbind(design, .padding = 0)
    B <- matrix(0, nrow = length(nuisance), ncol = length(l2))
    if (any(apply(design, 2L, stats::var) > 0) && stats::var(ystar) > 0) {
      path <- glmnet::glmnet(design, ystar, alpha = 1, intercept = FALSE,
                             standardize = TRUE, family = "gaussian")
      B[,] <- as.matrix(stats::coef(path, s = l2))[1L + seq_along(nuisance), , drop = FALSE]
    }
    partial <- matrix(a$Uy, nrow = length(a$d), ncol = length(l2)) - a$UZ %*% B
    out <- list(G1 = a$V %*% ((a$d / den) * partial), nuisance = B)
    if (keep) {
      out$debiased <- drop(a$V %*% (partial / a$d))
      out$eta <- sweep(a$U, 2L, a$d, "/") %*% t(a$V)
      if (nuisance.inference) {
        r <- a$y - drop(a$U %*% (a$d * crossprod(a$V, out$debiased)))
        out$ystar2 <- r - drop(a$U %*% (shrink * crossprod(a$U, r)))
        out$Zstar <- Zstar
      }
    }
    out
  }

  cv <- NULL
  if (tune) {
    error <- matrix(0, length(lambda1), length(lambda2))
    for (k in seq_len(max(foldid))) {
      valid <- which(foldid == k)
      a <- prepare(which(foldid != k))
      for (i in seq_along(lambda1)) {
        f <- fit.path(a, lambda1[i], lambda2)
        pred <- X[valid, G1, drop = FALSE] %*% f$G1 +
          X[valid, nuisance, drop = FALSE] %*% f$nuisance
        error[i, ] <- error[i, ] + colMeans((Y[valid] - pred)^2)
      }
    }
    cv <- expand.grid(lambda1 = lambda1, lambda2 = lambda2)
    cv$MSPE <- as.vector(error) / max(foldid)
    best <- which.min(cv$MSPE)
    lambda1 <- cv$lambda1[best]
    lambda2 <- cv$lambda2[best]
  }
  f <- fit.path(prepare(seq_len(n)), lambda1, lambda2, keep = TRUE)
  b <- stats::setNames(numeric(p), nm)
  b[G1] <- f$G1
  b[nuisance] <- f$nuisance
  bd <- b
  bd[G1] <- f$debiased
  intercept <- ymean - sum(xmean * b)
  intercept.d <- ymean - sum(xmean * bd)
  if (estimate.sigma) {
    variance.fit <- glmnet::cv.glmnet(X, Y, foldid = foldid, alpha = 1)
    slopes <- as.numeric(stats::coef(variance.fit, s = "lambda.min"))[-1L]
    df <- n - sum(slopes != 0)
    if (df <= 0) stop("Variance estimate has no residual degrees of freedom. Supply sigma.")
    sigma <- sqrt(sum((Y - X %*% slopes)^2) / df)
    if (!positive(sigma)) stop("Degenerate variance estimate. Supply sigma.")
  }
  tab <- NULL
  CI <- testpower <- NULL
  if (target.tests) {
    se <- sigma * sqrt(colSums(f$eta^2))
    z <- bd[G1] / se
    pval <- 2 * stats::pnorm(-abs(z))
    padj <- stats::p.adjust(pval, method = adjustment)
    level <- if (adjustment == "bonferroni") alpha / q else alpha
    width <- stats::qnorm(level / 2, lower.tail = FALSE) * se
    tab <- data.frame(index = G1, variable = nm[G1], estimate = b[G1],
                      debiased = bd[G1], se = se, lower = bd[G1] - width,
                      upper = bd[G1] + width, z = z, p.value = pval,
                      p.adjusted = padj, significant = padj <= alpha)
    CI <- as.matrix(tab[, c("lower", "upper"), drop = FALSE])
    testpower <- stats::setNames(as.integer(tab$significant), nm[G1])
  }
  nuisance.result <- NULL
  if (nuisance.inference) {
    if (!requireNamespace("hdi", quietly = TRUE)) stop("Please install hdi for nuisance inference.")
    if (is.null(nuisance.sigma)) nuisance.sigma <- sigma
    if (length(nuisance) < 2L) stop("The optional hdi step requires multiple nuisance columns.")
    nuisance.result <- hdi::lasso.proj(f$Zstar, f$ystar2, betainit = b[nuisance],
      sigma = nuisance.sigma, suppress.grouptesting = TRUE, verbose = FALSE)
  }
  list(call = match.call(), grouping = if (learned) "GCL" else "given",
       groups = stats::setNames(groups, nm), target = target, G1 = G1, G2d = nuisance,
       GCL.scores = scores, lambda1 = lambda1, lambda2 = lambda2,
       cv = cv, foldid = foldid, betaG1_hat = b[G1], betaG2d_hat = b[nuisance],
       betaG1_debiased = bd[G1], coefficients = b, intercept = intercept,
       coefficients.debiased = bd, intercept.debiased = intercept.d,
       fitted.values = ymean + drop(X %*% b), sigma = sigma,
       inference.method = if (target.tests) "Wald" else "none",
       alpha = alpha, adjustment = adjustment, inference.G1 = tab,
       nuisance.sigma = nuisance.sigma, nuisance.inference = nuisance.result,
       CI = CI, testpower = testpower)
}
