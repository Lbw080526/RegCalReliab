#' Naive Logistic Regression Analysis (Internal)
#'
#' Fits a logistic regression of the binary outcome \code{Y} on the subject-level
#' averaged exposure(s) \code{zbar}, with optional covariates \code{W.std}.
#' This baseline internal analysis ignores measurement-error correction beyond
#' taking replicate averages.
#'
#' @param zbar Numeric vector or matrix of standardized averages of exposure replicates.
#' @param W.std Optional standardized covariate matrix; default \code{NULL}.
#' @param Y Binary outcome vector of length n_m.
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
#'   \item Odds ratios (exponentiated coefficients)
#' }
#'
#' @noRd
#' @export
#' @importFrom stats glm vcov



naive_analysis_in_log = function(Y, zbar, W.std = NULL, sdz, sdw){


  if(is.null(W.std)){

  fit1 = glm(Y~zbar,family = "binomial")
  beta.fit1 = fit1$coefficients #naive estimator
  var1 = vcov(fit1) #naive covariance matrix

  tab1 = summary(fit1)$coefficients
  tab1[,1:2] = tab1[,1:2]/c(1,sdz)
  CI.low = tab1[,1]-1.96*tab1[,2]
  CI.high = tab1[,1]+1.96*tab1[,2]
  tab1 = cbind(tab1,exp(cbind(OR = tab1[, 1],CI.low,CI.high)))
  rownames(tab1) = sub("^zbar", "", rownames(tab1))

  }

  else{

    fit1 = glm(Y~ zbar + W.std,family = "binomial")
    beta.fit1 = fit1$coefficients #naive estimator
    var1 = vcov(fit1) #naive covariance matrix
    tab1 = summary(fit1)$coefficients
    tab1[,1:2] = tab1[,1:2]/c(1,sdz,sdw)
    CI.low = tab1[,1]-1.96*tab1[,2]
    CI.high = tab1[,1]+1.96*tab1[,2]
    tab1 = cbind(tab1,exp(cbind(OR = tab1[, 1],CI.low,CI.high)))
    rownames(tab1) = sub("^zbar", "", rownames(tab1))
    rownames(tab1) = sub("^W\\.std", "", rownames(tab1))
  }



  list(
    var1 = var1,
    `Naive estimates` = tab1
  )
}

