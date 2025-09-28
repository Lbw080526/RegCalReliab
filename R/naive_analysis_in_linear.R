#' Naive Linear Regression (Internal Reliability Study)
#'
#' \code{naive_analysis_in_linear()} fits a standard linear regression model
#' using subject-level averages of replicate measurements from an internal
#' reliability study—ignoring any further measurement-error correction.
#' It returns the uncorrected (naive) coefficient estimates, their standard
#' errors, confidence intervals, and the associated covariance matrix. These
#' results provide a benchmark for comparison with corrected regression
#' calibration estimates.
#'
#' @param zbar Numeric vector or matrix of standardized subject-level averages
#'   of exposure replicates (\eqn{n \times t}).
#' @param W.std Optional numeric matrix of standardized error-free covariates
#'   (\eqn{n \times q}); if not provided, the model is fit with exposures only.
#' @param Y Numeric outcome vector of length \eqn{n}.
#' @param sdz Numeric vector of length \eqn{t}, giving the standard deviations
#'   of the unstandardized exposures. Used to rescale coefficients back to the
#'   original measurement scale.
#' @param sdw Optional numeric vector of length \eqn{q}, giving the standard
#'   deviations of the unstandardized covariates. Used to rescale coefficients
#'   back to the original scale when \code{W.std} is supplied.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{var1}}{Covariance matrix of the naive linear regression estimates.}
#'   \item{\code{Naive estimates}}{Matrix of naive linear regression results,
#'         including coefficient estimates, standard errors, t-values,
#'         p-values, and 95\% confidence intervals on the original scale.}
#' }
#'
#' @details
#' This function is typically used as the first step in an internal reliability
#' study to illustrate the bias introduced by ignoring measurement error.
#' It scales coefficients and their standard errors back to the original scale
#' using the supplied standard deviations.
#'
#' @examples
#' set.seed(1)
#' # Simulated internal data: 80 subjects, 2 replicates of 1 exposure
#' z.rep <- cbind(rnorm(80), rnorm(80))
#' zbar <- rowMeans(z.rep)
#' Y <- 2 + 0.5 * zbar + rnorm(80)
#'
#' # Standardize and get SDs
#' zbar.std <- scale(zbar)
#' sdz <- sd(zbar)
#'
#' # Run naive linear regression ignoring measurement error
#' res <- naive_analysis_in_linear(
#'   zbar = zbar.std,
#'   W.std = NULL,
#'   Y = Y,
#'   sdz = sdz,
#'   sdw = NULL
#' )
#' str(res)
#'
#' @noRd



naive_analysis_in_linear = function(Y, zbar, W.std = NULL, sdz, sdw){

  z_df = as.data.frame(zbar)
  colnames(z_df) = colnames(zbar)

  if(is.null(W.std)){
  model_df = data.frame(Y = Y, z_df)
  fit1 = lm(Y ~ ., data = model_df)
  beta.fit1 = fit1$coefficients #naive estimator
  var1 = vcov(fit1) #naive covariance matrix

  tab1 = summary(fit1)$coefficients
  tab1[,1:2] = tab1[,1:2]/c(1,sdz)
  CI.low = tab1[,1]-1.96*tab1[,2]
  CI.high = tab1[,1]+1.96*tab1[,2]
  tab1[,1:2] <- tab1[,1:2] / c(1, sdz)

  }

  else{
    W_df = as.data.frame(W.std)
    colnames(W_df) = colnames(W.std)
    model_df = data.frame(Y = Y, z_df, W_df)
    fit1 = lm(Y ~ ., data = model_df)
    beta.fit1 = fit1$coefficients #naive estimator
    var1 = vcov(fit1) #naive covariance matrix
    tab1 = summary(fit1)$coefficients
    tab1[,1:2] = tab1[,1:2]/c(1,sdz,sdw)
    CI.low = tab1[,1]-1.96*tab1[,2]
    CI.high = tab1[,1]+1.96*tab1[,2]
    tab1[,1:2] <- tab1[,1:2] / c(1, sdz, sdw)
  }

  CI.low  <- tab1[,1] - 1.96 * tab1[,2]
  CI.high <- tab1[,1] + 1.96 * tab1[,2]
  tab1    <- cbind(tab1, CI.low = CI.low, CI.high = CI.high)



  list(
    var1 = var1,
    `Naive estimates` = tab1
  )
}

