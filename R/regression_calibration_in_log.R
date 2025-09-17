#' Regression Calibration for Logistic Regression (Internal Replicate Study)
#'
#' Performs regression calibration for a logistic GLM of the binary outcome \code{Y}
#' using subject-level replicate data from an internal reliability study. Corrects
#' for classical additive measurement error by estimating between-person and within-person
#' covariance components and computing calibrated predictors \code{xhat}.
#'
#' @param Y Binary outcome vector of length n.
#' @param zbar Matrix (n x t) of standardized subject-level averages of replicate exposures.
#' @param z.std List of length \code{t}; each element is an (n x r_i) matrix of
#'   standardized replicates for one exposure, padded with \code{NA} if fewer than
#'   the maximum number of replicates.
#' @param W.std Optional standardized covariate matrix (n x q); default \code{NULL}.
#' @param muz Mean(s) of exposures used for standardization.
#' @param muw Optional mean(s) of covariates used for standardization.
#' @param sdz Vector of standard deviations used to standardize exposures.
#' @param sdw Optional vector of standard deviations used to standardize covariates.
#' @param r Integer vector of replicate counts for each subject (length n).
#' @param var1 Naive covariance matrix of coefficients (from naive logistic model),
#'   used to set dimension names for covariance matrices.
#'
#' @return A list with:
#' \item{Corrected estimates}{Coefficient table with regression-calibrated estimates,
#'   robust (sandwich) standard errors, 95\% confidence intervals, and odds ratios.}
#' \item{icc}{Intraclass correlation matrix, quantifying reliability of the exposure(s).}
#' \item{sigmax}{Estimated between-person covariance matrix.}
#' \item{sigmawithin}{Estimated within-person covariance matrix.}
#' \item{sigmazstar}{Estimated total covariance matrix based on internal replicate structure.}
#' \item{sigmazhat}{Subject-specific total covariance matrices adjusted by replicate counts.}
#' \item{xhat}{Matrix of regression-calibrated predictors for the main study.}
#' \item{beta.fit2}{Vector of fitted coefficients from the corrected logistic regression.}
#' \item{v12star}{Between-person exposure block of the calibration matrix.}
#' \item{sigma}{Within-person variance of exposures estimated from replicates.}
#' \item{fit2}{Fitted \code{glm} object from the corrected logistic regression.}
#' \item{v}{Effective sample size adjustment factor for correlated replicates.}
#'
#' @details
#' Internal replicate data are used to decompose the variance of observed exposures
#' into between-person and within-person components. These components are combined
#' to compute regression-calibrated predictors \code{xhat}, which replace the
#' error-prone replicate averages \code{zbar} in the outcome model. Sandwich
#' variance estimators are used to obtain robust standard errors. When covariates
#' are included, their covariance with exposures is also incorporated into the
#' calibration matrix.
#'
#' @noRd
#' @export
#' @importFrom stats glm
#' @importFrom sandwich sandwich



reg_calibration_in_log = function(Y, zbar, z.std, W.std = NULL, muz, muw, sdz, sdw ,r, var1){
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

    fit2 = glm(Y~xhat,family = "binomial")
    beta.fit2 = fit2$coefficients
    var2 = sandwich::sandwich(fit2)

    tab2 = summary(fit2)$coefficients
    tab2[,2] = sqrt(diag(var2))
    tab2[,1:2] = tab2[,1:2]/c(1,sdz)
    CI.low = tab2[,1]-1.96*tab2[,2]
    CI.high = tab2[,1]+1.96*tab2[,2]
    tab2 = cbind(tab2,exp(cbind(OR = tab2[, 1],CI.low,CI.high)))
    rownames(tab2) = sub("^xhat", "", rownames(tab2))

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


    fit2 = glm(Y~xhat+W.std,family = "binomial")
    beta.fit2 = fit2$coefficients
    var2 = sandwich::sandwich(fit2)

    tab2 = summary(fit2)$coefficients
    tab2[,2] = sqrt(diag(var2))
    tab2[,1:2] = tab2[,1:2]/c(1,sdz,sdw)
    CI.low = tab2[,1]-1.96*tab2[,2]
    CI.high = tab2[,1]+1.96*tab2[,2]
    tab2 = cbind(tab2,exp(cbind(OR = tab2[, 1],CI.low,CI.high)))
    # Rename rows: clean variable names
    rownames(tab2) = sub("^xhat", "", rownames(tab2))
    rownames(tab2) <- sub("^W\\.std", "", rownames(tab2))


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



