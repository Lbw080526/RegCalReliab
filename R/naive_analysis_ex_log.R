#' Naive Logistic Regression Analysis (External)
#'
#' This function performs a logistic regression analysis without any measurement error correction.
#'
#' @param z.main.std Standardized main-study data.
#' @param W.main.std Standardized covariates for the main study (optional).
#' @param Y Binary outcome vector for the main-study subjects.
#' @param sdz Standard deviation of the exposures in the main study.
#' @param sdw Standard deviations of the covariates (optional).
#'
#' @return A list containing:
#' \itemize{
#'   \item{\code{var1}}: Covariance matrix for the naive logistic regression.
#'   \item{\code{Naive estimates}}: The naive logistic regression estimates and their confidence intervals.
#' }
#' @details This function runs a logistic regression using the main-study data without measurement error correction.
#' It returns both the estimates and the covariance matrix from the naive regression.
#' @noRd
#' @export

naive_analysis_ex_log = function(z.main.std, W.main.std = NULL, Y, sdz, sdw) {


  z_df = as.data.frame(z.main.std)
  colnames(z_df) = colnames(z.main.std)

  if(is.null(W.main.std)){

    model_df = data.frame(Y = Y, z_df)
    fit1 = glm(Y ~ ., data = model_df, family = "binomial")
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
    fit1 = glm(Y ~ ., data = model_df, family = "binomial") # Fit logistic regression using the error-prone exposure zmain and confounders Wmain
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
