#' Regression Calibration for Logistic Regression (External Reliability Study)
#'
#' This function applies regression calibration to correct measurement error in the logistic regression model.
#'
#' @param z.main.std Standardized main-study data.
#' @param z.rep.std Standardized reliability-study data.
#' @param r Integer vector of replicate counts for reliability-study subjects.
#' @param W.main.std Standardized covariates for the main study (optional).
#' @param Y Binary outcome vector for the main-study subjects.
#' @param muz Mean(s) of exposures used for standardization.
#' @param muw Optional mean(s) of covariates used for standardization.
#' @param sdz Standard deviations used to standardize z.main.
#' @param sdw Standard deviations used to standardize W.main.
#' @param indicator Indicator vector for main-study and reliability-study subjects.
#'
#' @return A list containing the regression calibration results, including calibrated estimates and variance components.
#' @noRd
#' @export




reg_calibration_ex_log = function(z.main.std, z.rep.std, r, W.main.std = NULL, Y, muz, muw, sdz, sdw, indicator) {

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
    fit2 = glm(Y~xhat,family = "binomial")
    beta.fit2 = fit2$coefficients
    var2 = sandwich::sandwich(fit2) # sandwich variance estimator

    tab2 = summary(fit2)$coefficients
    tab2[,2] = sqrt(diag(var2))
    tab2[,1:2] = tab2[,1:2]/c(1,sdz)
    CI.low = tab2[,1]-1.96*tab2[,2]
    CI.high = tab2[,1]+1.96*tab2[,2]
    tab2 = cbind(tab2,exp(cbind(OR = tab2[, 1],CI.low,CI.high)))

    rownames(tab2) <- sub("^xhat", "", rownames(tab2))


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
      fit2 = fit2,
      zbar = zbar
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
    fit2 = glm(Y~xhat + W.main.std,family = "binomial")
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
      fit2 = fit2,
      zbar = zbar
    ))

  }
}

