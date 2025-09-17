#' Naive Poisson (Log-Linear) Analysis (Internal)
#'
#' Fits a Poisson GLM with log link of the count outcome \code{Y} on the
#' subject-level averaged exposure(s) \code{zbar}, with optional covariates
#' \code{W.std}, without any measurement-error correction.
#'
#' @param zbar Numeric vector or matrix of standardized averages of exposure replicates.
#' @param W.std Optional standardized covariate matrix; default \code{NULL}.
#' @param Y Non-negative integer count outcome vector of length n_m.
#' @param sdz Vector of standard deviations used to standardize \code{zbar}.
#' @param sdw Optional vector of standard deviations used to standardize \code{W.std}.
#'
#'
#' @details
#' This baseline model ignores measurement error beyond averaging replicates
#' and is used for comparison against regression-calibrated analyses. The coefficient
#' table will include:
#' \itemize{
#'   \item Point estimates for coefficients
#'   \item Standard errors (SE)
#'   \item 95% Confidence intervals (CI) for each coefficient
#'   \item Exponentiated coefficients (rate ratios, RR)
#' }
#'
#' @noRd
#' @export
#' @importFrom stats glm vcov poisson


naive_analysis_in_poisson = function(Y, zbar, W.std = NULL, sdz, sdw) {

  if(is.null(W.std)) {
    # Fit Poisson regression with log link
    fit1 = glm(Y ~ zbar, family = poisson(link = "log"))
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
    rownames(tab1) = sub("^zbar", "", rownames(tab1))

  } else {
    # Fit Poisson regression with log link and covariates
    fit1 = glm(Y ~ zbar + W.std, family = poisson(link = "log"))
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
    rownames(tab1) = sub("^zbar", "", rownames(tab1))
    rownames(tab1) = sub("^W\\.std", "", rownames(tab1))
  }

  list(
    var1 = var1,
    `Naive estimates` = tab1
  )
}
