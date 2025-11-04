#' Naive Poisson Regression (External Reliability Study)
#'
#' \code{naive_analysis_ex_poisson()} fits a standard poisson regression model
#' using only the main-study data—ignoring any measurement error in the
#' exposure variables.
#' It returns the uncorrected (naive) coefficient estimates, their standard
#' errors, confidence intervals, exponentiated rate ratios, and the associated
#' covariance matrix. These results provide a benchmark for comparison with
#' corrected regression calibration estimates.
#'
#' @param z.main.std Numeric matrix of standardized main-study exposures
#'   (\eqn{n_m \times t}), typically the \code{z.main.std} output from
#'   \code{\link{prepare_data_ex}}.
#' @param W.main.std Optional numeric matrix of standardized error-free
#'   covariates (\eqn{n_m \times q}); if not provided, the model is fit with
#'   exposures only.
#' @param Y Non-negative integer count outcome vector of length \eqn{n_m}.
#' @param sdz Numeric vector of length \eqn{t}, giving the standard deviations
#'   of the unstandardized exposures. Used to rescale coefficients back to the
#'   original measurement scale.
#' @param sdw Optional numeric vector of length \eqn{q}, giving the standard
#'   deviations of the unstandardized covariates. Used to rescale coefficients
#'   back to the original scale when \code{W.main.std} is supplied.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{var1}}{Covariance matrix of the naive Poisson regression estimates.}
#'   \item{\code{Naive estimates}}{Matrix of naive Poisson regression results,
#'         including coefficient estimates, standard errors, z-values,
#'         p-values, exponentiated coefficients (rate ratios), and 95\% confidence
#'         intervals on the original scale.}
#' }
#'
#' @details
#' This function is typically used as the first step in an external reliability
#' study to illustrate the bias introduced by ignoring measurement error.
#' It scales coefficients and their standard errors back to the original scale
#' using the supplied standard deviations.
#'
#' @examples
#' set.seed(1)
#' # Simulated main-study data: 100 subjects, 1 exposure
#' z = matrix(rnorm(100), ncol = 1)
#' colnames(z) = "exposure"
#' Y = rpois(100, lambda = exp(0.2 + 0.3 * z))
#' sdz = apply(z, 2, sd)
#'
#' # Run naive Poisson regression ignoring measurement error
#' res = naive_analysis_ex_poisson(
#'   z.main.std = scale(z),
#'   W.main.std = NULL,
#'   Y = Y,
#'   sdz = sdz,
#'   sdw = NULL
#' )
#' str(res)
#'
#' @noRd


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
