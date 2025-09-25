#' Naive Poisson Regression (Internal Reliability Study)
#'
#' \code{naive_analysis_in_poisson()} fits a standard Poisson log-linear model
#' using only subject-level averaged replicates, ignoring measurement error
#' beyond replicate averaging. This provides the uncorrected (naive) estimates
#' as a benchmark for comparison with regression-calibrated results.
#'
#' @param zbar Numeric vector or matrix of standardized subject-level averaged
#'   exposures (\eqn{n \times t}), typically output from
#'   \code{\link{prepare_data_in}}.
#' @param W.std Optional numeric matrix of standardized error-free covariates
#'   (\eqn{n \times q}); default \code{NULL}.
#' @param Y Non-negative integer outcome vector of length \eqn{n}.
#' @param sdz Numeric vector of length \eqn{t}, giving the standard deviations
#'   of the unstandardized exposures (used to rescale coefficients).
#' @param sdw Optional numeric vector of length \eqn{q}, giving the standard
#'   deviations of the unstandardized covariates (used to rescale coefficients).
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{var1}}{Covariance matrix of the naive Poisson regression estimates.}
#'   \item{\code{Naive estimates}}{Matrix of naive Poisson regression results,
#'         including coefficient estimates, standard errors, z-values, p-values,
#'         rate ratios (RR), and 95\% confidence intervals on the original scale.}
#' }
#'
#' @details
#' This baseline model uses replicate-averaged exposures and ignores
#' measurement-error correction. Coefficients and standard errors are scaled
#' back to the original exposure (and covariate) scales using the supplied
#' standard deviations.
#'
#' @examples
#' set.seed(123)
#' # Simulated replicate data: 100 subjects, 1 exposure with 2 replicates
#' z.rep <- cbind(rnorm(100), rnorm(100))
#' zbar <- rowMeans(z.rep)
#' Y <- rpois(100, exp(0.3 + 0.5 * zbar))
#' sdz <- sd(zbar)
#'
#' # Run naive Poisson regression
#' res <- naive_analysis_in_poisson(
#'   Y = Y,
#'   zbar = scale(zbar),
#'   W.std = NULL,
#'   sdz = sdz,
#'   sdw = NULL
#' )
#' str(res)
#'
#' @noRd
#' @export



naive_analysis_in_poisson = function(Y, zbar, W.std = NULL, sdz, sdw) {

  z_df = as.data.frame(zbar)
  colnames(z_df) = colnames(zbar)

  if(is.null(W.std)) {
    # Fit Poisson regression with log link
    model_df = data.frame(Y = Y, z_df)
    fit1 = glm(Y ~ ., data = model_df, family = poisson(link = "log"))
    beta.fit1 = fit1$coefficients
    var1 = vcov(fit1)

    tab1 = summary(fit1)$coefficients
    # Scale coefficients and SEs
    tab1[,1:2] = tab1[,1:2]/c(1, sdz)
    # Calculate confidence intervals
    CI.low = tab1[,1] - 1.96*tab1[,2]
    CI.high = tab1[,1] + 1.96*tab1[,2]
    # Exponentiate to get rate ratios
    tab1 = cbind(tab1, exp(cbind(RR = tab1[, 1], CI.low, CI.high)))


  } else {
    # Fit Poisson regression with log link and covariates
    W_df = as.data.frame(W.std)
    colnames(W_df) = colnames(W.std)
    model_df = data.frame(Y = Y, z_df, W_df)
    fit1 = glm(Y ~ ., data = model_df, family = poisson(link = "log"))
    beta.fit1 = fit1$coefficients
    var1 = vcov(fit1)

    tab1 = summary(fit1)$coefficients
    # Scale coefficients and SEs
    tab1[,1:2] = tab1[,1:2]/c(1, sdz, sdw)
    # Calculate confidence intervals
    CI.low = tab1[,1] - 1.96*tab1[,2]
    CI.high = tab1[,1] + 1.96*tab1[,2]
    # Exponentiate to get rate ratios
    tab1 = cbind(tab1, exp(cbind(RR = tab1[, 1], CI.low, CI.high)))

  }

  list(
    var1 = var1,
    `Naive estimates` = tab1
  )
}
