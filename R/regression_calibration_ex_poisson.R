#' Regression Calibration for Poisson Regression (External Reliability Study)
#'
#' Performs regression calibration for a Poisson GLM with log link of the count outcome \code{Y}
#' using main-study error-prone exposures \code{z.main.std} and replicate data
#' from an external reliability study (\code{z.rep.std}, \code{r}, \code{indicator}).
#' Corrects for classical additive measurement error by estimating between-person
#' and within-person covariance components and computing calibrated predictors
#' \code{xhat}. Fits the corrected Poisson regression with optional covariates
#' \code{W.main.std}.
#'
#' @param z.main.std Matrix (n_m x t) of standardized error-prone exposure(s) from the main study.
#' @param z.rep.std List of length \code{t}; each element is an (n_r x r_i) matrix of
#'   standardized replicate measurements for the corresponding exposure in the reliability study.
#' @param r Integer vector of replicate counts for each subject (length n_m + n_r).
#'   Main-study subjects should have replicate count 1.
#' @param W.main.std Optional matrix (n_m x q) of standardized covariates for the main study; default \code{NULL}.
#' @param Y Non-negative integer count outcome vector of length n_m (main-study subjects).
#' @param muz Mean(s) of exposures used for standardization.
#' @param muw Optional mean(s) of covariates used for standardization.
#' @param sdz Vector of standard deviations used to standardize \code{z.main.std}.
#' @param sdw Optional vector of standard deviations used to standardize \code{W.main.std}.
#' @param indicator Binary vector (length n_m + n_r): 1 = main-study subject, 0 = reliability subject.
#'
#' @return A list with:
#' \item{Corrected estimates}{Coefficient table with regression-calibrated estimates,
#'   robust (sandwich) standard errors, 95\% confidence intervals, and exponentiated
#'   coefficients (rate ratios).}
#' \item{icc}{Intraclass correlation matrix, quantifying reliability of the exposure(s).}
#' \item{sigmax}{Estimated between-person covariance matrix.}
#' \item{sigmawithin}{Estimated within-person covariance matrix.}
#' \item{sigmaz}{Total observed covariance matrix of error-prone exposures (and covariates if included).}
#' \item{xhat}{Matrix of regression-calibrated predictors for the main study.}
#' \item{beta.fit2}{Vector of fitted coefficients from the corrected Poisson regression.}
#' \item{v12star}{Between-person exposure block of the calibration matrix.}
#' \item{sigma}{Within-person variance of exposures estimated from replicates.}
#' \item{fit2}{Fitted \code{glm} object from the corrected Poisson regression.}
#'
#' @details
#' The method uses replicate data to decompose the variance of observed exposures
#' into between-person and within-person components. These variance components
#' are then used to compute regression-calibrated predictors \code{xhat}, which
#' replace the error-prone exposures in the outcome model. Sandwich variance
#' estimators are used to obtain robust standard errors. When covariates are
#' included, their covariance with exposures is also incorporated into the
#' calibration matrix. Exponentiated coefficients are interpreted as rate ratios.
#'
#' @noRd
#' @export
#' @importFrom stats glm poisson
#' @importFrom sandwich sandwich



reg_calibration_ex_poisson = function(z.main.std, z.rep.std, r, W.main.std = NULL, Y, muz, muw, sdz, sdw, indicator) {

  # -----------------------------------------------
  # 0) Basic dimensions
  # -----------------------------------------------
  nm = nrow(z.main.std) # number of subjects in the main study
  nr = sum(indicator == 0) # number of subjects in the reliability study
  n = nm+nr # total number of subjects (main + reliability)
  t = ncol(z.main.std) # number of variables in z.main.std


  # -----------------------------------------------
  # 1) CASE 1:  W.std == NULL
  # -----------------------------------------------
  if (is.null(W.main.std)) {

    # Compute means for each subject from the reliability study (with replicates).
    zbar = sapply(z.rep.std,function(y) rowMeans(y,na.rm = T))

    r0 = r[indicator==0] # reliability study weights

    dif = sapply(1:nm, function(x) z.main.std[x,])

    # Calculate the covariance among zi, (error-prone meaurenment)
    if(t == 1){ # if just one variable in zmain , n*1, else n*m
      sigmaz = t(dif)%*%(dif)/nm #1*n * n*1
    }else{
      sigmaz = dif%*%t(dif)/nm #n*m * m*n
    }

    # From the reliability study, calculate the within-person variation.
    diff = as.matrix(na.omit(sapply(1:t,function(x) (z.rep.std[[x]]-zbar[,x]))))

    if(t==1){
      sigma = sum(sapply(1:nrow(diff),function(x) diff[x,]%*%t(diff[x,])))/sum(r0-1)
    }else{
      sigma = matrix(rowSums(sapply(1:nrow(diff),function(x) diff[x,]%*%t(diff[x,]))),t)/sum(r0-1)
    }

    # The matrix use to calibrate Z in X_hat,
    # extracts the between-person "exposure" block = (SigZ - Sigma_within)
    v12star = sigmaz[1:t,]-sigma

    # Compute the calibrated predictor, xhat. (For each main study subject,
    # we “correct” the observed measurement using the estimated variance components.)
    if (t == 1) {
      # Ensure matrix with 1 column
      xhat = matrix(apply(cbind(z.main.std), 1, function(tt)
        v12star %*% solve(sigmaz) %*% tt),
        ncol = 1)
    } else {
      xhat = t(apply(cbind(z.main.std), 1, function(tt)
        v12star %*% solve(sigmaz) %*% tt))
    }
    colnames(xhat) = colnames(z.main.std)



    # we get the between person variance (using sigamaz - sigma).
    sigmax = sigmaz-sigma
    sigmawithin = sigma
    icc = sigmax%*%solve(sigmaz)

    # -----  Fit Corrected Outcome Model -----
    fit2 = glm(Y ~ xhat, family = poisson(link = "log"))
    beta.fit2 = fit2$coefficients
    var2 = sandwich::sandwich(fit2) # sandwich variance estimator

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
      sigmaz = sigmaz,
      xhat = xhat,
      beta.fit2 = beta.fit2,
      v12star = v12star,
      sigma = sigma,
      fit2 = fit2
    ))

  }


  else {
    # -----------------------------------------------
    # 2) CASE 2:  W.std != NULL
    # -----------------------------------------------

    #restricted to the reliability subset
    nm = sum(indicator)
    nr = n-nm

    q = ncol(W.main.std)


    zbar = sapply(z.rep.std,function(y) rowMeans(y,na.rm = T))

    r0 = r[indicator==0]

    v = sum(r0)-sum(r0^2)/sum(r0) #scalar, measures how many independent data point when account for correlated replicates

    dif = rbind(sapply(1:nm, function(x) z.main.std[x,]),
                sapply(1:nm, function(x) W.main.std[x,]))

    sigmaz = dif%*%t(dif)/nm
    #matches \hat{\Sigma}_Z, \hat{\Sigma}_W, \hat{\Sigma}_{ZW}

    # matrix of "within-subject" differences from each replicate's subject-level mean
    diff = as.matrix(na.omit(sapply(1:t,function(x) (z.rep.std[[x]]-zbar[,x]))))

    # the estimated within-person variance
    if(t==1){
      sigma = sum(sapply(1:nrow(diff),function(x) diff[x,]%*%t(diff[x,])))/sum(r0-1)
    }else{
      sigma = matrix(rowSums(sapply(1:nrow(diff),function(x) diff[x,]%*%t(diff[x,]))),t)/sum(r0-1)
    }

    v12star = sigmaz[1:t,]-cbind(sigma,matrix(0,ncol = q,nrow=t))

    # Compute the corrected xhat
    if (t == 1) {
      xhat = matrix(apply(cbind(z.main.std, W.main.std), 1, function(tt)
        v12star %*% solve(sigmaz) %*% tt),
        ncol = 1)
    } else {
      xhat = t(apply(cbind(z.main.std, W.main.std), 1, function(tt)
        v12star %*% solve(sigmaz) %*% tt))
    }

    colnames(xhat) = colnames(z.main.std)
    colnames(W.main.std) = colnames(W.main.std)

    # -----  Fit Corrected Outcome Model -----

    # Fit logistic regression on the corrected exposures xhat + confounders Wmain
    fit2 = glm(Y ~ xhat + W.main.std, family = poisson(link = "log"))
    beta.fit2 = fit2$coefficients

    # A "sandwich" variance that partially adjusts for regression aspects but not fully for the measurement model
    var2 = sandwich::sandwich(fit2)

    # some parameters
    tab2 = summary(fit2)$coefficients
    tab2[,2] = sqrt(diag(var2))
    tab2[,1:2] = tab2[,1:2]/c(1,sdz,sdw)
    CI.low = tab2[,1]-1.96*tab2[,2]
    CI.high = tab2[,1]+1.96*tab2[,2]
    tab2 = cbind(tab2,exp(cbind(OR = tab2[, 1],CI.low,CI.high)))

    rownames(tab2) = sub("^xhat", "", rownames(tab2))

    rownames(tab2) = sub("^W.main.std", "", rownames(tab2))



    # sigmax = Sigma_X = total var - within-person var
    # entire between-person covariance matrix of dimension (t+q) * (t+q)
    sigmax = sigmaz-rbind(cbind(sigma,matrix(0,ncol = q,nrow=t)),matrix(0,ncol=t+q,nrow=q))
    sigmawithin = rbind(cbind(sigma,matrix(0,ncol = q,nrow=t)),matrix(0,ncol=t+q,nrow=q))

    # intraclass correlation matrix
    icc = sigmax%*%solve(sigmaz)

    return(list(
      `Corrected estimates` = tab2,
      icc = icc,
      sigmax = sigmax,
      sigmawithin = sigmawithin,
      sigmaz = sigmaz,
      xhat = xhat,
      beta.fit2 = beta.fit2,
      v12star = v12star,
      sigma = sigma,
      fit2 = fit2
    ))

  }
}

