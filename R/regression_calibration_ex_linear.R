#' Regression Calibration for Linear Regression (External Reliability Study)
#'
#' \code{reg_calibration_ex_linear()} corrects for classical additive measurement
#' error in logistic regression using data from an external reliability study.
#' It implements regression calibration by estimating between- and within-subject
#' variance components from replicate measurements and then fitting a logistic
#' regression to the calibrated exposures.
#' A robust (sandwich) variance estimator is used for valid inference.
#'
#' @param z.main.std Numeric matrix of standardized main-study exposures
#'   (\eqn{n_m \times t}), typically the \code{z.main.std} output from
#'   \code{\link{prepare_data_ex}}.
#' @param z.rep.std Named list of standardized replicate measurements from the
#'   external reliability study. Each list element is a matrix of dimension
#'   \eqn{n_r \times r_i} for a particular exposure.
#' @param r Integer vector of replicate counts for reliability-study subjects,
#'   of length \eqn{n_m + n_r}, with main-study subjects coded as 1.
#' @param W.main.std Optional numeric matrix of standardized error-free
#'   covariates (\eqn{n_m \times q}). If omitted, the calibration is performed
#'   for exposures only.
#' @param Y Numeric outcome vector of length \eqn{n_m}.
#' @param muz Numeric vector of means of the unstandardized exposures.
#' @param muw Optional numeric vector of means of the unstandardized covariates.
#' @param sdz Numeric vector of standard deviations of the unstandardized exposures.
#' @param sdw Optional numeric vector of standard deviations of the unstandardized covariates.
#' @param indicator Binary vector of length \eqn{n_m + n_r} indicating main-study
#'   (1) vs. reliability-study (0) subjects.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{Corrected estimates}}{Matrix of regression-calibrated
#'         coefficients, sandwich standard errors, t-values, p-values,
#'         and 95\% confidence intervals on the original scale.}
#'   \item{\code{icc}}{Intraclass correlation (matrix) quantifying reliability
#'         of the error-prone exposures.}
#'   \item{\code{sigmax}}{Estimated between-person covariance matrix of the true exposures.}
#'   \item{\code{sigmawithin}}{Estimated within-person (measurement-error) covariance matrix.}
#'   \item{\code{sigmaz}}{Estimated total covariance matrix of the observed exposures.}
#'   \item{\code{xhat}}{Matrix of calibrated exposure predictions used in the corrected regression.}
#'   \item{\code{beta.fit2}}{Vector of calibrated linear regression coefficients.}
#'   \item{\code{v12star}}{Calibration matrix used to map observed to corrected exposures.}
#'   \item{\code{sigma}}{Within-person variance matrix used in the calibration step.}
#'   \item{\code{fit2}}{The fitted \code{lm} object for the corrected linear regression.}
#' }
#'
#' @details
#' The method follows the classical regression calibration framework for
#' external reliability studies:
#' 1. Estimate total (\eqn{\Sigma_Z}) and within-subject (\eqn{\Sigma_\epsilon})
#'    covariance matrices using main-study and reliability-study data.
#' 2. Compute the between-subject covariance matrix (\eqn{\Sigma_X}) as
#'    \eqn{\Sigma_Z - \Sigma_\epsilon}.
#' 3. Calibrate each main-study measurement \eqn{Z_i} to
#'    \eqn{E[X_i | Z_i]} using the calibration matrix
#'    \eqn{\Sigma_X \Sigma_Z^{-1}}.
#' 4. Fit a linear regression of \eqn{Y} on \eqn{X^\text{hat}} (and optional covariates).
#'
#' @examples
#' set.seed(1)
#' # Simulated main-study data: 80 subjects, 1 exposure
#' z = matrix(rnorm(80), ncol = 1)
#' colnames(z) = "sbp"
#' Y = 2 + 0.5 * z + rnorm(80)
#'
#' # Reliability study: 40 subjects, 2 replicates
#' z.rep = list(sbp = matrix(rnorm(40 * 2), nrow = 40))
#' r = c(rep(1, 80), rep(2, 40))
#' indicator = c(rep(1, 80), rep(0, 40))
#'
#' # Standardize data
#' sdz = apply(z, 2, sd)
#' z.main.std = scale(z)
#' z.rep.std = list(sbp = scale(z.rep$sbp))
#'
#' # Run regression calibration
#' fit = reg_calibration_ex_linear(
#'   z.main.std = z.main.std,
#'   z.rep.std  = z.rep.std,
#'   r          = r,
#'   W.main.std = NULL,
#'   Y          = Y,
#'   muz        = colMeans(z),
#'   muw        = NULL,
#'   sdz        = sdz,
#'   sdw        = NULL,
#'   indicator  = indicator
#' )
#' str(fit)
#'
#' @noRd




reg_calibration_ex_linear = function(z.main.std, z.rep.std, r, W.main.std = NULL, Y, muz, muw, sdz, sdw, indicator) {

  # Basic dimensions
  nm = nrow(z.main.std) # number of subjects in the main study
  nr = sum(indicator == 0) # number of subjects in the reliability study
  n = nm+nr # total number of subjects (main + reliability)
  t = ncol(z.main.std) # number of variables in z.main.std


  # W.std == NULL
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
    if(t==1){
      xhat = apply(cbind(z.main.std),1,function(t) (v12star%*%solve(sigmaz)%*%t))
    }else{
      xhat = t(apply(cbind(z.main.std),1,function(t) (v12star%*%solve(sigmaz)%*%t)))
    }

    if (is.null(dim(xhat))) xhat = matrix(xhat, ncol = t)

    colnames(xhat) = colnames(z.main.std)

    # we get the between person variance (using sigamaz - sigma).
    sigmax = sigmaz-sigma
    sigmawithin = sigma
    icc = sigmax%*%solve(sigmaz)

    # Fit Corrected Outcome Model
    # Turn calibrated exposures into a data frame with correct names
    xhat_df = as.data.frame(xhat)
    colnames(xhat_df) = colnames(z.main.std)   # e.g. "sbp", "chol"
    model_df = data.frame(Y = Y, xhat_df)
    fit2 = lm(Y ~ ., data = model_df)
    beta.fit2 = fit2$coefficients
    var2 = sandwich::sandwich(fit2) # sandwich variance estimator

    tab2 = summary(fit2)$coefficients
    tab2[,2] = sqrt(diag(var2))
    tab2[,1:2] = tab2[,1:2]/c(1,sdz)
    CI.low = tab2[,1]-1.96*tab2[,2]
    CI.high = tab2[,1]+1.96*tab2[,2]
    tab2 = cbind(tab2, CI.low = CI.low, CI.high = CI.high)

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
    #W.std != NULL
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
    if(t==1){
      xhat = apply(cbind(z.main.std, W.main.std),1,function(t) (v12star%*%solve(sigmaz)%*%t))
    }else{
      xhat = t(apply(cbind(z.main.std, W.main.std),1,function(t) (v12star%*%solve(sigmaz)%*%t)))
    }

    if (is.null(dim(xhat))) xhat = matrix(xhat, ncol = t)  # <<< FIX
    colnames(xhat) = colnames(z.main.std)
    colnames(W.main.std) = colnames(W.main.std)

    # Fit Corrected Outcome Model

    xhat_df = as.data.frame(xhat)
    W_df = as.data.frame(W.main.std)
    colnames(xhat_df) = colnames(z.main.std)   # e.g. "sbp", "chol"
    colnames(W_df) = colnames(W.main.std)
    model_df = data.frame(Y = Y, xhat_df, W_df)
    fit2 = lm(Y ~ ., data = model_df)
    beta.fit2 = fit2$coefficients

    # A "sandwich" variance that partially adjusts for regression aspects but not fully for the measurement model
    var2 = sandwich::sandwich(fit2)

    # some parameters
    tab2 = summary(fit2)$coefficients
    tab2[,2] = sqrt(diag(var2))
    tab2[,1:2] = tab2[,1:2]/c(1,sdz,sdw)
    CI.low = tab2[,1]-1.96*tab2[,2]
    CI.high = tab2[,1]+1.96*tab2[,2]
    tab2 = cbind(tab2, CI.low = CI.low, CI.high = CI.high)

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

