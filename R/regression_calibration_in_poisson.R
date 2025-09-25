#' Regression Calibration for Poisson Regression (Internal Reliability Study)
#'
#' \code{reg_calibration_in_poisson()} corrects for classical additive
#' measurement error in Poisson log-linear regression using replicate
#' measurements from an internal reliability study. It estimates between- and
#' within-subject variance components from replicates, constructs calibrated
#' predictors \code{xhat}, and fits a Poisson GLM to the calibrated exposures.
#' A robust (sandwich) variance estimator is used for valid inference.
#'
#' @param Y Integer, non-negative count outcome vector of length \eqn{n}.
#' @param zbar Numeric matrix (\eqn{n \times t}) of standardized subject-level
#'   averages of replicate exposures.
#' @param z.std Named list of length \eqn{t}; each element is an \eqn{n \times r_i}
#'   matrix of standardized replicate measurements for one exposure, padded
#'   with \code{NA} if fewer than the maximum number of replicates.
#' @param W.std Optional numeric matrix of standardized error-free covariates
#'   (\eqn{n \times q}). If omitted, calibration is performed for exposures only.
#' @param muz Numeric vector of means of the unstandardized exposures.
#' @param muw Optional numeric vector of means of the unstandardized covariates.
#' @param sdz Numeric vector of standard deviations of the unstandardized exposures.
#' @param sdw Optional numeric vector of standard deviations of the unstandardized covariates.
#' @param r Integer vector of replicate counts for each subject (length \eqn{n}).
#' @param var1 Naive covariance matrix of coefficients (from
#'   \code{\link{naive_analysis_in_poisson}}), used to set dimension names for
#'   covariance matrices.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{Corrected estimates}}{Matrix of calibrated Poisson regression
#'         coefficients, robust (sandwich) standard errors, z-values, p-values,
#'         rate ratios (RR), and 95\% confidence intervals on the original scale.}
#'   \item{\code{icc}}{Intraclass correlation (matrix) quantifying reliability
#'         of the error-prone exposures.}
#'   \item{\code{sigmax}}{Estimated between-person covariance matrix of the true exposures.}
#'   \item{\code{sigmawithin}}{Estimated within-person (measurement-error) covariance matrix.}
#'   \item{\code{sigmazstar}}{Estimated total covariance matrix based on replicate structure.}
#'   \item{\code{sigmazhat}}{List of subject-specific total covariance matrices
#'         adjusted for replicate counts.}
#'   \item{\code{xhat}}{Matrix of calibrated exposure predictions used in the corrected regression.}
#'   \item{\code{beta.fit2}}{Vector of calibrated Poisson regression coefficients.}
#'   \item{\code{v12star}}{Calibration matrix used to map observed to corrected exposures.}
#'   \item{\code{sigma}}{Within-person variance matrix estimated from replicates.}
#'   \item{\code{fit2}}{The fitted \code{glm} object for the corrected Poisson regression.}
#'   \item{\code{v}}{Effective sample size adjustment factor for correlated replicates.}
#' }
#'
#' @details
#' The method follows the classical regression calibration framework for
#' internal reliability studies:
#' \enumerate{
#'   \item Estimate total (\eqn{\Sigma_Z}) and within-subject (\eqn{\Sigma_\epsilon})
#'         covariance matrices using replicate data.
#'   \item Compute the between-subject covariance matrix (\eqn{\Sigma_X}) as
#'         \eqn{\Sigma_Z - \Sigma_\epsilon}.
#'   \item Construct subject-specific total covariance matrices adjusted for
#'         replicate counts.
#'   \item Calibrate each subject’s replicate average \eqn{Z_i} to
#'         \eqn{E[X_i \mid Z_i]} using the calibration matrix.
#'   \item Fit a Poisson log-linear model using the calibrated exposures
#'         \eqn{X_i^\text{hat}} and compute sandwich standard errors.
#' }
#'
#' @examples
#' set.seed(123)
#' # Internal reliability study: 60 subjects, 2 replicates of 1 exposure
#' z.rep <- cbind(rnorm(60), rnorm(60))
#' zbar <- rowMeans(z.rep)
#' Y <- rpois(60, lambda = exp(0.2 + 0.5 * zbar))
#'
#' # Standardize data
#' zbar.std <- scale(zbar)
#' sdz <- sd(zbar)
#' z.std <- list(sbp = scale(z.rep))
#' r <- rep(2, 60)
#'
#' # Naive covariance (for dimension labels)
#' naive <- naive_analysis_in_poisson(Y = Y, zbar = zbar.std,
#'                                    W.std = NULL, sdz = sdz, sdw = NULL)
#'
#' # Apply regression calibration
#' fit <- reg_calibration_in_poisson(
#'   Y = Y,
#'   zbar = as.matrix(zbar.std),
#'   z.std = z.std,
#'   W.std = NULL,
#'   muz = mean(zbar),
#'   muw = NULL,
#'   sdz = sdz,
#'   sdw = NULL,
#'   r = r,
#'   var1 = naive$var1
#' )
#' str(fit)
#'
#' @noRd
#' @export





reg_calibration_in_poisson = function(Y, zbar, z.std, W.std = NULL, muz, muw, sdz, sdw ,r, var1){
  # -----------------------------------------------
  # 0) Basic dimensions
  # -----------------------------------------------
  n = length(r)
  t = length(z.std)
  q = ncol(W.std)

  # -----------------------------------------------
  # 1) CASE 1:  W.std == NULL
  # -----------------------------------------------
  if(is.null(W.std)){

    v = sum(r)-sum(r^2)/sum(r)

    dif = sapply(1:n, function(x) zbar[x,])

    if(t == 1){
      sigmazstar = t(dif)%*%(dif*r)/v
    }else{
      sigmazstar = dif%*%t(dif*r)/v
    }

    diff = as.matrix(na.omit(sapply(1:t,function(x) (z.std[[x]]-zbar[,x]))))
    if(t==1){
      sigma = sum(sapply(1:nrow(diff),function(x) diff[x,]%*%t(diff[x,])))/sum(r-1)
    }else{
      sigma = matrix(rowSums(sapply(1:nrow(diff),function(x) diff[x,]%*%t(diff[x,]))),t)/sum(r-1)
    }
    v12star = sigmazstar - (n)*sigma/v
    sigmax = sigmazstar-(n)*sigma/v
    sigmawithin = sigma
    icc = sigmax%*%solve(sigmazstar)

    if(t==1){
      sigmazhat = t(sapply(r,function(x) sigmazstar-(n)*sigma/v+sigma/x))
    }else{
      sigmazhat = sapply(r,function(x) sigmazstar-(n)*sigma/v+sigma/x)
    }
    if(t==1){
      xhat = sapply(1:n,function(i) (v12star%*%solve(matrix(sigmazhat[,i],ncol=t))%*%zbar[i]))
      xhat <- matrix(xhat, ncol = 1)                 # <- new
      colnames(xhat) <- colnames(zbar) %||% "z"
    }else{
      xhat = t(sapply(1:n,function(i) (v12star%*%solve(matrix(sigmazhat[,i],ncol=t))%*%zbar[i,])))
    }

    colnames(xhat) = paste0(colnames(zbar))

    xhat_df = as.data.frame(xhat)
    colnames(xhat_df) = colnames(zbar)   # e.g. "sbp", "chol"
    model_df = data.frame(Y = Y, xhat_df)
    fit2 = glm(Y ~ ., data = model_df, family = poisson(link = "log"))
    beta.fit2 = fit2$coefficients
    var2 = sandwich::sandwich(fit2)

    tab2 = summary(fit2)$coefficients
    tab2[,2] = sqrt(diag(var2))
    tab2[,1:2] = tab2[,1:2]/c(1,sdz)
    CI.low = tab2[,1]-1.96*tab2[,2]
    CI.high = tab2[,1]+1.96*tab2[,2]
    tab2 = cbind(tab2,exp(cbind(OR = tab2[, 1],CI.low,CI.high)))

    return(list(
      `Corrected estimates` = tab2,
      icc = icc,
      sigmax = sigmax,
      sigmawithin= sigmawithin,
      xhat = xhat,
      beta.fit2 = beta.fit2,
      v12star = v12star,
      sigma = sigma,
      fit2 = fit2,
      sigmazstar = sigmazstar,
      sigmazhat = sigmazhat,
      v = v

    ))


  }


  else {
    # -----------------------------------------------
    # 2) CASE 2:  W.std != NULL
    # -----------------------------------------------
    v = sum(r)-sum(r^2)/sum(r)

    dif = rbind(sapply(1:n, function(x) zbar[x,]),
                sapply(1:n, function(x) W.std[x,]))


    sigmazstar = dif%*%t(dif*r)/v

    diff = as.matrix(na.omit(sapply(1:t,function(x) (z.std[[x]]-zbar[,x]))))

    if(t==1){
      sigma = sum(sapply(1:nrow(diff),function(x) diff[x,]%*%t(diff[x,])))/sum(r-1)
    }else{
      sigma = matrix(rowSums(sapply(1:nrow(diff),function(x) diff[x,]%*%t(diff[x,]))),t)/sum(r-1)
    }
    v12star = sigmazstar[1:t,]-cbind((n)*sigma/v,matrix(0,ncol = q,nrow=t))

    sigmax = sigmazstar-rbind(cbind((n)*sigma/v,matrix(0,ncol = q,nrow=t)),matrix(0,ncol=t+q,nrow=q))
    sigmawithin = rbind(cbind(sigma,matrix(0,ncol = q,nrow=t)),matrix(0,ncol=t+q,nrow=q))
    icc = sigmax%*%solve(sigmazstar)
    colnames(sigmax) = colnames(var1)[-1]
    rownames(sigmax) = rownames(var1)[-1]
    colnames(sigmawithin) = colnames(var1)[-1]
    rownames(sigmawithin) = rownames(var1)[-1]
    colnames(icc) = colnames(var1)[-1]
    rownames(icc) = rownames(var1)[-1]

    sigmazhat = sapply(r,function(x) sigmax+rbind(cbind(sigma/x,matrix(0,ncol = q,nrow=t)),matrix(0,ncol=t+q,nrow=q)))

    if(t==1){
      xhat = sapply(1:n,function(i) (v12star%*%solve(matrix(sigmazhat[,i],ncol=t+q))%*%cbind(zbar,W.std)[i,]))
      xhat <- matrix(xhat, ncol = 1)                 # <- new
      colnames(xhat) <- colnames(zbar) %||% "z"
    }else{
      xhat = t(sapply(1:n,function(i) (v12star%*%solve(matrix(sigmazhat[,i],ncol=t+q))%*%cbind(zbar,W.std)[i,])))
    }

    colnames(xhat) = paste0(colnames(zbar))
    colnames(W.std) = colnames(W.std)

    xhat_df <- as.data.frame(xhat)
    W_df <- as.data.frame(W.std)
    colnames(xhat_df) <- colnames(xhat)   # e.g. "sbp", "chol"
    colnames(W_df) <- colnames(W.std)
    model_df <- data.frame(Y = Y, xhat_df, W_df)
    fit2 <- glm(Y ~ ., data = model_df, family = poisson(link = "log"))
    beta.fit2 = fit2$coefficients
    var2 = sandwich::sandwich(fit2)

    tab2 = summary(fit2)$coefficients
    tab2[,2] = sqrt(diag(var2))
    tab2[,1:2] = tab2[,1:2]/c(1,sdz,sdw)
    CI.low = tab2[,1]-1.96*tab2[,2]
    CI.high = tab2[,1]+1.96*tab2[,2]
    tab2 = cbind(tab2,exp(cbind(OR = tab2[, 1],CI.low,CI.high)))


    return(list(
      `Corrected estimates` = tab2,
      icc = icc,
      sigmax = sigmax,
      sigmawithin = sigmawithin,
      xhat = xhat,
      beta.fit2 = beta.fit2,
      v12star = v12star,
      sigma = sigma,
      fit2 = fit2,
      sigmazstar = sigmazstar,
      sigmazhat = sigmazhat,
      v = v
    ))

  }

}



