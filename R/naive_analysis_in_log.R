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

  z_df = as.data.frame(zbar)
  colnames(z_df) = colnames(zbar)

  if(is.null(W.std)){
  model_df = data.frame(Y = Y, z_df)
  fit1 = glm(Y ~ ., data = model_df, family = "binomial")
  beta.fit1 = fit1$coefficients #naive estimator
  var1 = vcov(fit1) #naive covariance matrix

  tab1 = summary(fit1)$coefficients
  tab1[,1:2] = tab1[,1:2]/c(1,sdz)
  CI.low = tab1[,1]-1.96*tab1[,2]
  CI.high = tab1[,1]+1.96*tab1[,2]
  tab1 = cbind(tab1,exp(cbind(OR = tab1[, 1],CI.low,CI.high)))

  }

  else{
    W_df = as.data.frame(W.std)
    colnames(W_df) = colnames(W.std)
    model_df = data.frame(Y = Y, z_df, W_df)
    fit1 = glm(Y ~ ., data = model_df, family = "binomial")
    beta.fit1 = fit1$coefficients #naive estimator
    var1 = vcov(fit1) #naive covariance matrix
    tab1 = summary(fit1)$coefficients
    tab1[,1:2] = tab1[,1:2]/c(1,sdz,sdw)
    CI.low = tab1[,1]-1.96*tab1[,2]
    CI.high = tab1[,1]+1.96*tab1[,2]
    tab1 = cbind(tab1,exp(cbind(OR = tab1[, 1],CI.low,CI.high)))
  }



  list(
    var1 = var1,
    `Naive estimates` = tab1
  )
}

