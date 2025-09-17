#' Naive Poisson Regression Analysis (External)
#'
#' Fits a Poisson GLM of the count outcome \code{Y} on the observed
#' (error-prone) exposure(s) \code{z.main.std}, with optional covariates
#' \code{W.main.std}, without any measurement-error correction.
#'
#' @param z.main.std Matrix (n_m x t) of standardized error-prone exposure(s) from the main study.
#' @param W.main.std Optional matrix (n_m x q) of standardized covariates from the main study; default \code{NULL}.
#' @param Y Non-negative integer count outcome vector of length n_m.
#' @param sdz Vector of standard deviations used to standardize \code{z.main.std}.
#' @param sdw Optional vector of standard deviations used to standardize \code{W.main.std}.
#'
#'
#' @details
#' This baseline model ignores measurement error and is used for comparison
#' against regression-calibrated analyses. The coefficient table includes:
#' \itemize{
#'   \item Point estimates for coefficients
#'   \item Standard errors (SE)
#'   \item 95\% confidence intervals (CI)
#'   \item Exponentiated coefficients (rate ratios)
#' }
#'
#' @noRd
#' @export
#' @importFrom stats glm vcov poisson


naive_analysis_ex_poisson = function(z.main.std, W.main.std = NULL, Y, sdz, sdw) {


  z_df = as.data.frame(z.main.std)
  colnames(z_df) = colnames(z.main.std)

  if(is.null(W.main.std)){

    model_df = data.frame(Y = Y, z_df)
    fit1 = glm(Y ~ ., data = model_df, family = poisson(link = "log"))
    beta.fit1 = (fit1$coefficients) #naive estimator
    var1 = vcov(fit1) #naive covariance matrix

    # Adjust coefficient and standard error scales (dividing by sdz)
    tab1 = summary(fit1)$coefficients
    tab1[,1:2] = tab1[,1:2]/c(1,sdz)
    CI.low = tab1[,1]-1.96*tab1[,2]
    CI.high = tab1[,1]+1.96*tab1[,2]
    tab1 = cbind(tab1,exp(cbind(OR = tab1[, 1],CI.low,CI.high)))
  }

  else{

    W_df = as.data.frame(W.main.std)
    colnames(W_df) = colnames(W.main.std)
    model_df = data.frame(Y = Y, z_df, W_df)
    fit1 = glm(Y ~ ., data = model_df, family = poisson(link = "log"))
    beta.fit1 = (fit1$coefficients) #naive estimator
    var1 = vcov(fit1) #naive covariance matrix

    tab1 = summary(fit1)$coefficients
    tab1[,1:2] = tab1[,1:2]/c(1,sdz,sdw)
    CI.low = tab1[,1]-1.96*tab1[,2]
    CI.high = tab1[,1]+1.96*tab1[,2]
    tab1 = cbind(tab1,exp(cbind(OR = tab1[, 1],CI.low,CI.high))) # Combine everything
  }

  list(
    var1 = var1,
    `Naive estimates` = tab1
  )

}
