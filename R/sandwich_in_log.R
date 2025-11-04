#' Sandwich Variance Estimator for Internal Logistic Regression Calibration
#'
#' \code{sandwich_estimator_in_log()} computes robust (sandwich) standard errors
#' for regression calibration in logistic regression using internal replicate
#' measurements. It propagates uncertainty from estimating calibration
#' parameters (between- and within-subject variance components) into the final
#' coefficient variances, yielding valid confidence intervals and p-values.
#'
#' @param xhat Numeric matrix of calibrated exposures (\eqn{n \times t}),
#'   typically from \code{\link{reg_calibration_in_log}}.
#' @param zbar Numeric matrix (\eqn{n \times t}) of standardized subject-level
#'   replicate averages.
#' @param z.std Named list of standardized replicate exposure matrices.
#'   Each element is an \eqn{n \times r_i} matrix, possibly padded with \code{NA}.
#' @param r Integer vector of replicate counts per subject (length \eqn{n}).
#' @param Y Binary (0/1) outcome vector of length \eqn{n}.
#' @param v12star Calibration submatrix (between-person exposure block).
#' @param beta.fit2 Numeric vector of logistic regression coefficients from the
#'   calibrated model (\code{fit2}).
#' @param W.std Optional numeric matrix of standardized error-free covariates
#'   (\eqn{n \times q}); if provided, covariate variability is incorporated
#'   into the sandwich correction.
#' @param sigma Within-person variance matrix of replicate exposures.
#' @param sigmawithin Estimated within-person covariance matrix.
#' @param sigmazstar Total covariance matrix of standardized exposures
#'   (and covariates, if present).
#' @param sigmazhat List (or array) of subject-specific estimated covariance
#'   matrices adjusted for replicate counts.
#' @param sdz Numeric vector of standard deviations of the unstandardized exposures.
#' @param sdw Optional numeric vector of standard deviations of the unstandardized covariates.
#' @param muz Numeric vector of means of the unstandardized exposures.
#' @param muw Optional numeric vector of means of the unstandardized covariates.
#' @param fit2 The fitted \code{glm} object returned by \code{reg_calibration_in_log()}.
#' @param v Normalization constant used in the variance scaling step.
#'
#' @return A list with one component:
#' \describe{
#'   \item{\code{Sandwich Corrected estimates}}{Matrix of regression calibration
#'         estimates with sandwich-corrected standard errors, z-values,
#'         p-values, odds ratios (OR), and 95\% confidence intervals on the
#'         original scale.}
#' }
#'
#' @details
#' This function implements the full sandwich variance correction for internal
#' reliability studies:
#' \enumerate{
#'   \item \strong{Variance components}: estimates within- and between-subject
#'         covariance from replicate data.
#'   \item \strong{Calibration}: computes calibrated predictors
#'         \eqn{X^\text{hat}} using the estimated covariance components.
#'   \item \strong{Outcome model}: fits a logistic regression and applies the
#'         sandwich (robust) variance formula to propagate uncertainty from
#'         both calibration and outcome modeling.
#' }
#'
#' @examples
#' set.seed(123)
#' # Simulate internal replicate data: 50 subjects, 2 replicates of 1 exposure
#' z.rep = cbind(rnorm(50), rnorm(50))
#' zbar = rowMeans(z.rep)
#' Y = rbinom(50, 1, plogis(0.5 * zbar))
#'
#' # Standardize
#' zbar.std = scale(zbar)
#' sdz = sd(zbar)
#' z.std = list(sbp = scale(z.rep))
#' r = rep(2, 50)
#'
#' # (In practice, run reg_calibration_in_log() first to obtain xhat, fit2, etc.)
#' # Here we show a simplified call with mock objects:
#' # sandwich_estimator_in_log(xhat, zbar.std, z.std, r, Y,
#' #   v12star, beta.fit2, W.std = NULL, sigma, sigmawithin,
#' #   sigmazstar, sigmazhat, sdz, NULL, muz, NULL, fit2, v)
#'
#' @references
#' Carroll RJ, Ruppert D, Stefanski LA, Crainiceanu CM. \emph{Measurement Error
#' in Nonlinear Models: A Modern Perspective}. Chapman & Hall/CRC, 2006.
#' White H. A heteroskedasticity-consistent covariance matrix estimator and a
#' direct test for heteroskedasticity. \emph{Econometrica}. 1980;48(4):817–838.
#'
#' @noRd





sandwich_estimator_in_log = function(xhat,zbar,z.std,r,Y,v12star,beta.fit2,W.std = NULL,sigma,
                                     sigmawithin,sigmazstar,sigmazhat, sdz,sdw,muz,muw,fit2,v){

  # Basic dimensions
  n = length(r)
  t = length(z.std)
  q = ncol(W.std)

  muz_std = colMeans(zbar, na.rm = TRUE)  # length t

  if (!is.null(W.std)) {
    muw_std = colMeans(W.std, na.rm = TRUE)  # length q
  } else {
    muw_std = NULL
  }


  if(is.null(W.std)){

    p = as.vector(exp(beta.fit2 %*% t(cbind(1,xhat)))/(1+exp(beta.fit2 %*% t(cbind(1,xhat)))))
    sigmaz.inv = solve(sigmazstar)

    Z = t(cbind(zbar))

    m = matrix(0,nrow = t,ncol = t)
    c = matrix(0,nrow = t,ncol = t)

    ### sigmax
    ###diag
    ddv.x = sapply(1:t,function(x){
      a = rep(0,t)
      a[x] = a[x]+1
      diag(a)
    },simplify = F)
    ddm.x = sapply(1:t,function(x){
      a = rep(0,t)
      a[x] = a[x]+1
      diag(a)
    },simplify = F)
    db.x = sapply(1:t,function(x) t(sapply(1:n,function(y) (ddv.x[[x]]%*%solve(matrix(sigmazhat[,y],ncol=t))-
                                                              v12star%*%solve(matrix(sigmazhat[,y],ncol=t))%*%ddm.x[[x]]%*%solve(matrix(sigmazhat[,y],ncol=t)))%*%Z[,y])),simplify = F)

    ###off-diag
    if(t>1){odv.x = sapply(1:(t-1),function(x) sapply(min((x+1),t):t,function(y){
      c[x,y] = c[x,y]+1
      c[y,x] = c[y,x]+1
      c
    },simplify = F),simplify = F)
    odm.x = sapply(1:(t-1),function(x) sapply(min((x+1),t):t,function(y){
      m[x,y] = m[x,y]+1
      m[y,x] = m[y,x]+1
      m
    },simplify = F),simplify = F)
    ob.x = sapply(1:(t-1),function(x) sapply(1:(t-x),function(y)
      t(sapply(1:n,function(u) (odv.x[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t))-
                                  v12star%*%solve(matrix(sigmazhat[,u],ncol=t))%*%odm.x[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t)))%*%Z[,u])),simplify = F),simplify = F)
    }

    ### sigma
    ###diag
    db.0 = sapply(1:t,function(x) t(sapply(1:n, function(u) t((-ddv.x[[x]]%*%solve(matrix(sigmazhat[,u],ncol=t)))%*%(Z/r)[,u]))),simplify = F)

    ###off-diag
    if(t>1){ob.0 = sapply(1:(t-1),function(x) sapply(1:(t-x),function(y)
      t(sapply(1:n,function(u) t((-odv.x[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t)))%*%(Z/r)[,u]))),simplify = F),simplify = F)}

    m = (t)*(t+1)/2+(t*(t+1))/2 #number of the covariance estimates
    s = t+1 #number of beta estimates

    b = if(t>1){
      list(db.x,ob.x,db.0,ob.0)
    }else list(db.x,db.0)

    b = sapply(1:m,function(x) matrix(unlist(b),nrow=n)[,(t*x-t+1):(t*x)],simplify = F)

    d = as.matrix(p*(1-p)%*%t(beta.fit2[2:(t+1)]))/n

    bstar = sapply(b,function(x) -rowSums(x*d))

    a = matrix(NA,ncol=m,nrow=s)
    a[1,] = colSums(bstar)
    a[2:(t+1),] = t(apply(as.matrix(xhat), 2, function(x) colSums(x * bstar)))

    w = diag(p * (1 - p))
    Astar = t(cbind(1,xhat)) %*% w %*% cbind(1,xhat) /n


    A = rbind(cbind(diag(rep(-1,m)),
                    matrix(0,nrow=m,ncol=s)),cbind(a,Astar))
    A = rbind(cbind(matrix(diag(rep(-1,t),nrow=t),t),matrix(0,nrow=t,ncol=m+s)),
              cbind(matrix(0,nrow=m+s,ncol=t),A))

    dmgi = sapply(1:(t),function(x) r*(((Z[x,])-muz_std[x]))/sum(r))
    ###sigmas
    ###diag
    dgi = r*t(((Z-muz_std)^2-diag(sigmazstar))/v)
    ###off-diag
    if(t>1){ogi = sapply(1:(t-1),function(x) sapply((x+1):(t), function(y)
      (r*((Z[x,]-muz_std[x])*((Z[y,]-muz_std[y])-sigmazstar[x,y]))/v)))}
    ###within variance
    ###diag
    dgi.0 = if(t==1){
      sapply(1:t,function(y) (rowSums((z.std[[y]]-zbar[,y])^2,na.rm = T)-sigma*(r-1))/sum(r-1))
    }else sapply(1:t,function(y) (rowSums((z.std[[y]]-zbar[,y])^2,na.rm = T)-diag(sigma)[y]*(r-1))/sum(r-1))
    ###off-diag
    if(t>1){ogi.0 = sapply(1:(t-1),function(x) sapply((x+1):t, function(y)
      (rowSums((z.std[[x]]-zbar[,x])*(z.std[[y]]-zbar[,y]),na.rm = T)-sigma[x,y]*(r-1))/sum(r-1)))
    }

    ###betas
    gi.beta = cbind(1,xhat)*(Y-p) / n


    gi.1 = if(t==1){
      cbind(dmgi,dgi,gi.beta)
    }else cbind(dmgi,dgi,matrix(unlist(ogi),nrow=n),gi.beta)
    B1 = t(gi.1)%*%gi.1

    gi.2 = if(t==1){
      dgi.0
    }else cbind(dgi.0,matrix(unlist(ogi.0),nrow=n))
    B2 = t(gi.2)%*%gi.2

    m1 = (t)*(t+1)/2
    m2 = t*(t+1)/2
    B = rbind(cbind(B1[1:(t+m1),1:(t+m1)],matrix(0,t+m1,m2),B1[1:(t+m1),(t+m1+1):(t+m1+s)]),
              cbind(matrix(0,m2,t+m1),B2,matrix(0,m2,s)),
              cbind(B1[(t+m1+1):(t+m1+s),1:(t+m1)],matrix(0,s,m2),B1[(t+m1+1):(t+m1+s),(t+m1+1):(t+m1+s)]))

    cov2 = solve(A)%*%B%*%t(solve(A))

    tab3 = summary(fit2)$coefficients
    tab3[,2] = sqrt(diag(cov2)[(t+m+1):(t+m+s)])
    tab3[,1:2] = tab3[,1:2]/c(1,sdz)
    CI.low = tab3[,1]-1.96*tab3[,2]
    CI.high = tab3[,1]+1.96*tab3[,2]
    tab3[,3] = tab3[,1]/tab3[,2]
    tab3[,4] = 2*pnorm(tab3[,3],lower.tail=FALSE)
    tab3 = cbind(tab3,exp(cbind(OR = tab3[, 1],CI.low,CI.high)))

    return(list(
      `Sandwich Corrected estimates` = tab3))



  } else{

    p = as.vector(exp(beta.fit2 %*% t(cbind(1,xhat,W.std)))/(1+exp(beta.fit2 %*% t(cbind(1,xhat,W.std)))))
    sigmaz.inv = solve(sigmazstar)

    Z = t(cbind(zbar,W.std))

    m = matrix(0,nrow = t+q,ncol = t+q)
    c = matrix(0,nrow = t,ncol = t+q)

    ### sigmax
    ###diag
    ddv.x = sapply(1:t,function(x){
      a = rep(0,t)
      a[x] = a[x]+1
      cbind(diag(a),matrix(0,nrow = t,ncol = q))
    },simplify = F)
    ddm.x = sapply(1:t,function(x){
      a = rep(0,t+q)
      a[x] = a[x]+1
      diag(a)
    },simplify = F)
    db.x = sapply(1:t,function(x) t(sapply(1:n,function(y) (ddv.x[[x]]%*%solve(matrix(sigmazhat[,y],ncol=t+q))-
                                                              v12star%*%solve(matrix(sigmazhat[,y],ncol=t+q))%*%ddm.x[[x]]%*%solve(matrix(sigmazhat[,y],ncol=t+q)))%*%Z[,y])),simplify = F)

    ###off-diag
    if(t>1){odv.x=sapply(1:(t-1),function(x) sapply(min((x+1),t):t,function(y){
      c[x,y] = c[x,y]+1
      c[y,x] = c[y,x]+1
      c
    },simplify = F),simplify = F)
    odm.x = sapply(1:(t-1),function(x) sapply(min((x+1),t):t,function(y){
      m[x,y] = m[x,y]+1
      m[y,x] = m[y,x]+1
      m
    },simplify = F),simplify = F)

    ob.x = sapply(1:(t-1),function(x) sapply(1:(t-x),function(y)
      t(sapply(1:n,function(u) (odv.x[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t+q))-
                                  v12star%*%solve(matrix(sigmazhat[,u],ncol=t+q))%*%odm.x[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t+q)))%*%Z[,u])),simplify = F),simplify = F)
    }

    ### sigmaxW
    ### all off-diag
    odv.xw = sapply(1:t,function(x) sapply(max((x+1),t+1):(t+q),function(y){
      c[x,y] = c[x,y]+1
      c
    },simplify = F),simplify = F)
    odm.xw = sapply(1:t,function(x) sapply(max((x+1),t+1):(t+q),function(y){
      m[x,y] = m[x,y]+1
      m[y,x] = m[y,x]+1
      m
    },simplify = F),simplify = F)

    ob.xw = sapply(1:t,function(x) sapply(1:q,function(y)
      t(sapply(1:n,function(u) t((odv.xw[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t+q))-
                                    v12star%*%solve(matrix(sigmazhat[,u],ncol=t+q))%*%odm.xw[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t+q)))%*%Z[,u]))),simplify = F),simplify = F)

    ### sigmaW
    ###diag
    ddm.w = sapply((t+1):(t+q),function(x){
      a = rep(0,t+q)
      a[x] = a[x]+1
      diag(a)
    },simplify = F)
    db.w = sapply(1:q,function(x)
      t(sapply(1:n,function(u) t((-v12star%*%solve(matrix(sigmazhat[,u],ncol=t+q))%*%ddm.w[[x]]%*%solve(matrix(sigmazhat[,u],ncol=t+q)))%*%Z[,u]))),simplify = F)

    ###off-diag
    if(q>1){odm.w=sapply((t+1):(t+q-1),function(x) sapply(min((x+1),t+q):(t+q),function(y){
      m[x,y] = m[x,y]+1
      m[y,x] = m[y,x]+1
      m
    },simplify = F),simplify = F)
    ob.w = sapply(1:(q-1),function(x) sapply(1:(q-x),function(y)
      t(sapply(1:n,function(u) (-v12star%*%solve(matrix(sigmazhat[,u],ncol=t+q))%*%odm.w[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t+q)))%*%Z[,u])),simplify = F),simplify = F)
    }
    ### sigma
    ###diag
    db.0 = sapply(1:t,function(x) t(sapply(1:n, function(u) t((-ddv.x[[x]]%*%solve(matrix(sigmazhat[,u],ncol=t+q)))%*%(Z/r)[,u]))),simplify = F)

    ###off-diag
    if(t>1){ob.0 = sapply(1:(t-1),function(x) sapply(1:(t-x),function(y)
      t(sapply(1:n,function(u) t((-odv.x[[x]][[y]]%*%solve(matrix(sigmazhat[,u],ncol=t+q)))%*%(Z/r)[,u]))),simplify = F),simplify = F)}

    m = (t+q)*(t+q+1)/2+(t*(t+1))/2 #number of the covariance estimates
    s = t+q+1 #number of beta estimates

    b = if(t>1 & q>1){
      list(db.x,db.w,ob.x,ob.xw,ob.w,db.0,ob.0)
    }else if(q>1){
      list(db.x,db.w,ob.xw,ob.w,db.0)
    }else if(t>1){
      list(db.x,db.w,ob.x,ob.xw,db.0,ob.0)
    }else list(db.x,db.w,ob.xw,db.0)

    b = sapply(1:m,function(x) matrix(unlist(b),nrow=n)[,(t*x-t+1):(t*x)],simplify = F)

    d = as.matrix(p*(1-p)%*%t(beta.fit2[2:(t+1)]))/n

    bstar = sapply(b,function(x) -rowSums(x*d))

    a = matrix(NA,ncol=m,nrow=s)
    a[1,] = colSums(bstar)
    a[2:(t+1),] = t(apply(as.matrix(xhat), 2, function(x) colSums(x * bstar)))
    a[(t+2):(t+q+1),] = t(apply(as.matrix(W.std), 2, function(x) colSums(x * bstar)))

    w = diag(p * (1 - p))
    Astar = t(cbind(1,xhat, W.std)) %*% w %*% cbind(1,xhat, W.std)/n


    A = rbind(cbind(diag(rep(-1,m)),
                    matrix(0,nrow=m,ncol=s)),cbind(a,Astar))
    A = rbind(cbind(diag(rep(-1,t+q),nrow=t+q),matrix(0,nrow=t+q,ncol=m+s)),
              cbind(matrix(0,nrow=m+s,ncol=t+q),A))

    ###sigmas
    ###diag
    mu  = c(muz_std, muw_std)
    dmgi = sapply(1:(t+q),function(x) r*(Z[x,]-mu[x])/sum(r))

    dgi = sapply(1:(t+q),function(x) (r*((Z[x,]-mu[x])^2-diag(sigmazstar)[x]))/v)
    ###off-diag
    ogi = sapply(1:(t+q-1),function(x) sapply((x+1):(t+q), function(y)
      (r*((Z[x,]-mu[x])*(Z[y,]-mu[y])-sigmazstar[x,y]))/v))

    ###within variance
    ###diag
    dgi.0 = if(t==1){
      sapply(1:t,function(y) (rowSums((z.std[[y]]-zbar[,y])^2,na.rm = T)-sigma*(r-1))/sum(r-1))
    }else sapply(1:t,function(y) (rowSums((z.std[[y]]-zbar[,y])^2,na.rm = T)-diag(sigma)[y]*(r-1))/sum(r-1))
    ###off-diag
    if(t>1){ogi.0 = sapply(1:(t-1),function(x) sapply((x+1):t, function(y)
      (rowSums((z.std[[x]]-zbar[,x])*(z.std[[y]]-zbar[,y]),na.rm = T)-sigma[x,y]*(r-1))/sum(r-1)))
    }

    ###betas
    gi.beta = cbind(1,xhat,W.std)*(Y-p)/n

    gi.1 = cbind(dmgi,dgi,matrix(unlist(ogi),nrow=n),gi.beta)
    B1 = crossprod(gi.1)

    gi.2 = if(t==1){
      dgi.0
    }else cbind(dgi.0,matrix(unlist(ogi.0),nrow=n))
    B2 = t(gi.2)%*%gi.2

    m1 = (t+q)*(t+q+1)/2
    m2 = t*(t+1)/2
    B = rbind(cbind(B1[1:(t+q+m1),1:(t+q+m1)],matrix(0,t+q+m1,m2),B1[1:(t+q+m1),(t+q+m1+1):(t+q+m1+s)]),
              cbind(matrix(0,m2,t+q+m1),B2,matrix(0,m2,s)),
              cbind(B1[(t+q+m1+1):(t+q+m1+s),1:(t+q+m1)],matrix(0,s,m2),B1[(t+q+m1+1):(t+q+m1+s),(t+q+m1+1):(t+q+m1+s)]))

    cov2 = solve(A)%*%B%*%t(solve(A))

    CI.low = beta.fit2-1.96*sqrt(diag(cov2)[(t+q+m+1):(t+q+m+s)])
    CI.high = beta.fit2+1.96*sqrt(diag(cov2)[(t+q+m+1):(t+q+m+s)])
    coverage.xhat0 = log(1.5)<=CI.high & log(1.5)>=CI.low

    tab3 = summary(fit2)$coefficients
    tab3[,2] = sqrt(diag(cov2)[(t+q+m+1):(t+q+m+s)])
    tab3[,1:2] = tab3[,1:2]/c(1,sdz,sdw)
    CI.low = tab3[,1]-1.96*tab3[,2]
    CI.high = tab3[,1]+1.96*tab3[,2]
    tab3[,3] = tab3[,1]/tab3[,2]
    tab3[,4] = 2*pnorm(tab3[,3],lower.tail=FALSE)
    tab3 = cbind(tab3,exp(cbind(OR = tab3[, 1],CI.low,CI.high)))

    return(list(
      `Sandwich Corrected estimates` = tab3))


  }


}
