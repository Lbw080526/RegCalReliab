#' Sandwich Variance Estimator for External Poisson (Log-Linear) Regression
#'
#' Computes robust (sandwich) variance-corrected coefficient estimates for
#' Poisson log-linear regression models using external reliability study data.
#' Accounts for measurement error by incorporating both main-study and
#' replicate reliability study information into calibration and variance
#' estimation.
#'
#' @param xhat Matrix of regression-calibrated predictors for the main study.
#' @param z.main.std Standardized main-study exposure matrix (n_main x t).
#' @param z.rep.std List of standardized replicate exposure matrices for the
#'   reliability study.
#' @param r Integer vector of replicate counts for each subject.
#' @param Y Integer, non-negative count outcome vector.
#' @param indicator Binary indicator vector: \code{1 = main study}, \code{0 = reliability study}.
#' @param v12star Calibration submatrix (between-person exposure block).
#' @param beta.fit2 Estimated regression coefficients from the corrected Poisson model.
#' @param W.main.std Optional standardized covariate matrix (n_main x q); default \code{NULL}.
#' @param sigma Within-person covariance matrix of replicate exposures.
#' @param sigmaz Total covariance matrix of exposures (and covariates if present).
#' @param sigmawithin Estimated within-person covariance matrix.
#' @param sdz Vector of standard deviations used for exposure standardization.
#' @param sdw Optional vector of standard deviations used for covariate standardization.
#' @param muz Mean(s) of exposures used in standardization.
#' @param muw Optional mean(s) of covariates used in standardization.
#' @param fit2 Fitted Poisson regression model object from the calibration step.
#'
#' @return A list containing:
#' \item{Sandwich Corrected estimates}{Coefficient table with regression-calibrated
#'   estimates, robust (sandwich) standard errors, Wald test statistics, 95% confidence
#'   intervals, and exponentiated rate ratios (RR).}
#'
#' @details
#' The estimator augments Poisson regression calibration with external replicate data.
#' It partitions variance into between- and within-subject components, and builds
#' the full sandwich covariance estimator by combining contributions from main and
#' reliability studies. Derivative blocks for covariance parameters are constructed
#' to form the Jacobian and score variance matrices. This ensures valid inference
#' under classical additive measurement error when external reliability data are
#' available.
#'
#' @noRd
#' @export
#' @importFrom stats glm pnorm poisson



sandwich_estimator_ex_poisson = function(xhat,z.main.std,z.rep.std,r,Y,indicator, v12star,beta.fit2,W.main.std = NULL,sigma,
                                         sigmaz,sigmawithin,sdz,sdw,muz,muw,fit2){



  # -----------------------------------------------
  # 0) Basic dimensions
  # -----------------------------------------------
  nm = nrow(z.main.std) # number of subjects in the main study
  nr = sum(indicator == 0)
  n = nm+nr # total number of subjects (main + reliability)
  t = ncol(z.main.std) # number of variables in z.main.std

  muz_std <- as.numeric(colMeans(z.main.std))            # length t
  if (!is.null(W.main.std)) muw_std <- as.numeric(colMeans(W.main.std))  # length q


  if(is.null(W.main.std)){

    # Compute fitted probabilities p from the corrected (calibrated) logistic model
    p = as.vector(exp(beta.fit2 %*% t(cbind(1, xhat))))

    # Compute the inverse of the estimated covariance matrix of Z, sigmaz
    sigmaz.inv = solve(sigmaz)

    # Create matrix Z as the transpose of the main study data, from nm x t to t x nm
    Z = t(cbind(z.main.std))

    # Initialize two t x t zero matrices, they will be used later
    # to construct off–diagonal derivative matrices.
    m = matrix(0,nrow = t,ncol = t)
    c = matrix(0,nrow = t,ncol = t)
    ##sigmaX
    ###diag

    # ddv.x and ddm.x are essentially indicator matrices with a 1 in the (x,x) entry. S and T
    ddv.x = sapply(1:t,function(x){
      a = rep(0,t)
      a[x] = a[x]+1
      diag(a)
    },simplify = FALSE)

    ddm.x = sapply(1:t,function(x){
      a = rep(0,t)
      a[x] = a[x]+1
      diag(a)
    },simplify = FALSE)

    # the matrix of partial derivatives (with respect to the (x,x) diagonal element of sigmaz,) of the calibration
    # Diagonal Partials wrt sigmaz
    db.x = sapply(1:t,function(x) t((ddv.x[[x]]%*%sigmaz.inv-
                                       v12star%*%sigmaz.inv%*%ddm.x[[x]]%*%sigmaz.inv)%*%Z),simplify = FALSE)

    ###off-diag
    # Build off–diagonal derivative matrices
    # For each x from 1 to t-1, for each y from x+1 to t, update matrix c

    if(t>1){
      odv.x<-sapply(1:(t-1),function(x) sapply(min((x+1),t):t,function(y){
        c[x,y]<-c[x,y]+1
        c[y,x]<-c[y,x]+1
        c
      },simplify = FALSE),simplify = FALSE)
      odm.x<-sapply(1:(t-1),function(x) sapply(min((x+1),t):t,function(y){
        m[x,y]<-m[x,y]+1
        m[y,x]<-m[y,x]+1
        m
      },simplify = FALSE),simplify = F)

      # Compute the off-diagonal derivative contributions
      ob.x<-sapply(1:(t-1),function(x) sapply(1:(t-x),function(y)
        t((odv.x[[x]][[y]]%*%sigmaz.inv-
             v12star%*%sigmaz.inv%*%odm.x[[x]][[y]]%*%sigmaz.inv)%*%Z),simplify = FALSE),simplify = FALSE)
    }

    # Off‐Diagonal Partials wrt sigmaz


    ##sigma (the within–group variance)
    #Partials wrt Another Variance W
    ###diag
    db.0 = sapply(1:t,function(x) t((ddv.x[[x]]%*%sigmaz.inv)%*%Z),simplify = FALSE)
    # correspond to partial derivatives wrt the diagonal of the other covariance matrix sigmaZW

    ###off-diag
    if(t>1){ob.0 = sapply(1:(t-1),function(x) sapply(1:(t-x),function(y)
      t((odv.x[[x]][[y]]%*%sigmaz.inv)%*%Z),simplify = FALSE),simplify = FALSE)}
    # correspond to partial derivatives wrt the off-diagonal of the other covariance matrix sigmaZW


    m = (t)*(t+1)/2+(t*(t+1))/2 #number of the covariance estimates
    s = t+1 #number of beta estimates (regression parameters)

    # If t > 1, include both diagonal and off-diagonal pieces otherwise, only the diagonal pieces
    # When t=1 we have only one “error‐prone” covariate—so there are no off‐diagonal elements
    # (1×1 covariance matrix has no off‐diagonal entries)
    b = if(t>1){
      list(db.x,ob.x,db.0,ob.0)
    }else list(db.x,db.0)

    # Reshape the list b into a list of m matrices
    # For each index from 1 to m, extract a matrix with nm rows and t columns from the unlisted b
    b = sapply(1:m,function(x) matrix(unlist(b),nrow=nm)[,(t*x-t+1):(t*x)],simplify = FALSE)

    #d = as.matrix(p*(1-p)%*%t(beta.fit2[2:(t+1)])) / nm
    d = as.matrix(p %*% t(beta.fit2[2:(t+1)])) / nm  # CHANGED: Removed (1-p)

    bstar = sapply(b,function(x) -rowSums(x*d)) # - sum_j b_{ij} * d_{ij}

    a = matrix(NA,ncol=m,nrow=s)
    a[1,] = colSums(bstar)
    a[2:(t+1),] = t(apply(as.matrix(xhat), 2, function(x) colSums(x * bstar))) #correspond to A2


    w = diag(p)  # Poisson variance is p

    Astar = t(cbind(1,xhat)) %*% w %*% cbind(1,xhat) /nm # correspond to A3

    A = rbind(cbind(diag(rep(-1,m)),
                    matrix(0,nrow=m,ncol=s)),cbind(a,Astar)) # building up A2 and A3 together

    A = rbind(cbind(matrix(diag(rep(-1,t),nrow=t),t),matrix(0,nrow=t,ncol=m+s)),
              cbind(matrix(0,nrow=m+s,ncol=t),A)) # combining A1, A2, A3 together building the A matrix

    ###sigmas
    ###diag
    dmgi = sapply(1:t, function(x) (Z[x,] - muz_std[x]) / nm)

    #dgi<-sapply(1:(t),function(x) (((Z[x,]-mu[x])^2-diag(sigmaz)[x]))/nm)

    dgi = t(((Z - muz_std)^2 - diag(sigmaz)) / nm)

    ###off-diag
    if (t > 1) {
      ogi <- sapply(1:(t-1), function(x) sapply((x+1):t, function(y)
        ((Z[x,] - muz_std[x]) * (Z[y,] - muz_std[y]) - sigmaz[x,y]) / nm))
    }

    ###within variance
    ###diag
    zbar = sapply(z.rep.std,function(y) rowMeans(y,na.rm = T))
    r0 = r[indicator==0]

    dgi.0<-if(t==1){
      sapply(1:t,function(y) (rowSums((z.rep.std[[y]]-zbar[,y])^2,na.rm = T)-sigma*(r0-1))/sum(r0-1))
    }else sapply(1:t,function(y) (rowSums((z.rep.std[[y]]-zbar[,y])^2,na.rm = T)-diag(sigma)[y]*(r0-1))/sum(r0-1))

    ###off-diag
    if(t>1){ogi.0<-sapply(1:(t-1),function(x) sapply((x+1):t, function(y)
      (rowSums((z.rep.std[[x]]-zbar[,x])*(z.rep.std[[y]]-zbar[,y]),na.rm = T)-sigma[x,y]*(r0-1))/sum(r0-1)))
    }


    ###betas
    gi.beta = cbind(1,xhat)*(Y-p) / nm # Construct the residualfrom the logistic model.


    # Combine the residuals for the main‐study portion into one matrix
    gi.1 = if(t==1){
      cbind(dmgi,dgi,gi.beta)
    }else cbind(dmgi,dgi,matrix(unlist(ogi),nrow=nm),gi.beta) # adding the off-diagonal

    B1 = t(gi.1)%*%gi.1 #the sigmaV_iV_i^T pattern but only for the main study

    gi.2 = if(t==1){
      dgi.0
    }else cbind(dgi.0,matrix(unlist(ogi.0),nrow=nr))# Combine the residuals for the reliability study into one matrix gi.2

    B2 = t(gi.2)%*%gi.2 # the sigmaV_iV_i^T pattern but only for the reliable study

    m1 = (t)*(t+1)/2 # number of unique covariance parameters (lower‐triangular elements) in main study
    m2 = t*(t+1)/2 # number of unique covariance parameters (lower‐triangular elements) in reliable study

    # Block [1 : (t+m1), 1 : (t+m1)], from B1 (main study)
    # Block (2,2): B2 = reliability part.
    # Block [ (t+m1+1):(t+m1+s), (t+m1+1):(t+m1+s) ] from B1 = logistic part.
    B = rbind(cbind(B1[1:(t+m1),1:(t+m1)],matrix(0,t+m1,m2),B1[1:(t+m1),(t+m1+1):(t+m1+s)]),
              cbind(matrix(0,m2,t+m1),B2,matrix(0,m2,s)),
              cbind(B1[(t+m1+1):(t+m1+s),1:(t+m1)],matrix(0,s,m2),B1[(t+m1+1):(t+m1+s),(t+m1+1):(t+m1+s)]))

    cov2 = solve(A)%*%B%*%t(solve(A)) # Sandwich Computation
    var3 = cov2[(t+m+1):(t+m+s),(t+m+1):(t+m+s)] # Extract the submatrix corresponding to the logistic coefficients

    tab3 = summary(fit2)$coefficients
    tab3[,2] = sqrt(diag(cov2)[(t+m+1):(t+m+s)])
    tab3[,1:2] = tab3[,1:2]/c(1,sdz)
    CI.low = tab3[,1]-1.96*tab3[,2]
    CI.high = tab3[,1]+1.96*tab3[,2]
    tab3[,3] = tab3[,1]/tab3[,2]  #  Z‐value
    tab3[,4] = 2*pnorm(tab3[,3],lower.tail=FALSE) # two‐sided p‐value from the normal distribution
    tab3 = cbind(tab3,exp(cbind(OR = tab3[, 1],CI.low,CI.high)))

    rownames(tab3) = sub("^xhat", "", rownames(tab3))

    return(list(
      `Sandwich Corrected estimates` = tab3
    ))


  } else{

    colnames(xhat) = colnames(z.main.std)
    colnames(W.main.std) = colnames(W.main.std)
    nm = sum(indicator)
    nr = n-nm
    q = ncol(W.main.std)

    model_matrix = cbind(1, xhat, W.main.std)
    colnames(model_matrix) = names(beta.fit2)
    #p = as.vector(exp(beta.fit2 %*% t(model_matrix))/(1+exp(beta.fit2 %*% t(model_matrix))))
    p = as.vector(exp(beta.fit2 %*% t(model_matrix)))

    # invert the total covariance matrix of (zmain, Wmain)
    # the inv. in the passage, use to correct for the dependencies between Z and W when
    # computing sensitivity measures. It appear in the sandwich esitmator to ensure proper variance
    sigmaz.inv = solve(sigmaz)

    # the transpose of the design matrix for the original measured
    #  lines up with computing sum(Zi-mu_z) or similar terms in the score equations
    Z = t(cbind(z.main.std,W.main.std))

    # initialize zero matrices used to mark partial-derivative positions for off-diagonal blocks.
    m = matrix(0,nrow = t+q,ncol = t+q)
    c = matrix(0,nrow = t,ncol = t+q)

    ##sigmaX
    ###diag

    # S matrix for diagonal elements of the top-left block sigmaz
    ddv.x = sapply(1:t,function(x){
      a = rep(0,t)
      a[x] = a[x]+1
      cbind(diag(a),matrix(0,nrow = t,ncol = q))
    },simplify = FALSE)

    # T matrix for diagonal elements
    ddm.x = sapply(1:t,function(x){
      a = rep(0,t+q)
      a[x]<-a[x]+1
      diag(a)
    },simplify = FALSE)

    # compute the final partial derivative derivaritve x_hat w.r.t to sigmaz using the formula
    # bi = (S*Inv - v*Inv*T*Inv)*Ci
    db.x = sapply(1:t,function(x) t((ddv.x[[x]]%*%sigmaz.inv-
                                       v12star%*%sigmaz.inv%*%ddm.x[[x]]%*%sigmaz.inv)%*%Z),simplify = FALSE)

    # Off-Diagonal version:
    ###off-diag
    if(t>1){odv.x = sapply(1:(t-1),function(x) sapply(min((x+1),t):t,function(y){
      c[x,y] = c[x,y]+1
      c[y,x] = c[y,x]+1
      c
    },simplify = FALSE),simplify = FALSE)
    odm.x = sapply(1:(t-1),function(x) sapply(min((x+1),t):t,function(y){
      m[x,y] = m[x,y]+1
      m[y,x] = m[y,x]+1
      m
    },simplify = FALSE),simplify = FALSE)
    ob.x = sapply(1:(t-1),function(x) sapply(1:(t-x),function(y)
      t((odv.x[[x]][[y]]%*%sigmaz.inv-
           v12star%*%sigmaz.inv%*%odm.x[[x]][[y]]%*%sigmaz.inv)%*%Z),simplify = FALSE),simplify = FALSE)
    }

    ##sigmaXW (ZW)
    ##all off-diag
    # S for the cross block
    odv.xw = sapply(1:t,function(x) sapply(max((x+1),t+1):(t+q),function(y){
      c[x,y]<-c[x,y]+1
      c
    },simplify = FALSE),simplify = FALSE)
    # T
    odm.xw = sapply(1:t,function(x) sapply(max((x+1),t+1):(t+q),function(y){
      m[x,y]<-m[x,y]+1
      m[y,x]<-m[y,x]+1
      m
    },simplify = FALSE),simplify = FALSE)

    # Cross-Covariance, sigma ZW
    ob.xw = sapply(1:t,function(x) sapply(1:q,function(y)
      t((odv.xw[[x]][[y]]%*%sigmaz.inv-
           v12star%*%sigmaz.inv%*%odm.xw[[x]][[y]]%*%sigmaz.inv)%*%Z),simplify = FALSE),simplify = FALSE)

    ##sigmaW
    ###diag
    ddm.w = sapply((t+1):(t+q),function(x){
      a = rep(0,t+q)
      a[x] = a[x]+1
      diag(a)
    },simplify = FALSE)

    # sigmaW
    db.w = sapply(1:q,function(x)
      t((-v12star%*%sigmaz.inv%*%ddm.w[[x]]%*%sigmaz.inv)%*%Z),simplify = FALSE)

    ###off-diag
    if(q>1){odm.w = sapply((t+1):(t+q-1),function(x) sapply(min((x+1),t+q):(t+q),function(y){
      m[x,y]<-m[x,y]+1
      m[y,x]<-m[y,x]+1
      m
    },simplify = FALSE),simplify = FALSE)

    # the bi formula
    ob.w = sapply(1:(q-1),function(x) sapply(1:(q-x),function(y)
      t((-v12star%*%sigmaz.inv%*%odm.w[[x]][[y]]%*%sigmaz.inv)%*%Z),simplify = FALSE),simplify = FALSE)
    }

    ##sigma within-person sigma
    ###diag
    db.0<-sapply(1:t,function(x) t((-ddv.x[[x]]%*%sigmaz.inv)%*%Z),simplify = FALSE)

    ###off-diag
    if(t>1){ob.0 = sapply(1:(t-1),function(x) sapply(1:(t-x),function(y)
      t((-odv.x[[x]][[y]]%*%sigmaz.inv)%*%Z),simplify = FALSE),simplify = FALSE)}

    m = (t+q)*(t+q+1)/2+(t*(t+1))/2 #number of the covariance estimates
    s = t+q+1 #number of beta estimates
    #put partial derivatives into a list b, base on t and q
    b = if(t>1 & q>1){
      list(db.x,db.w,ob.x,ob.xw,ob.w,db.0,ob.0)
    }else if(q>1){
      list(db.x,db.w,ob.xw,ob.w,db.0)
    }else if(t>1){
      list(db.x,db.w,ob.x,ob.xw,db.0,ob.0)
    }else list(db.x,db.w,ob.xw,db.0)
    #Reformat all derivative blocks into m submatrices, each corresponding to a specific covariance parameter.


    b = sapply(1:m,function(x) matrix(unlist(b),nrow=nm)[,(t*x-t+1):(t*x)],simplify = FALSE)

    # logistic derivative w.r.t. X: p_i(1-p_i) * betaZ
    # matches di = p_i(1-p_i)*betaZ in the passage
    #d = as.matrix(p*(1-p)%*%t(beta.fit2[2:(t+1)]))/nm
    d = as.matrix(p %*% t(beta.fit2[2:(t+1)])) / nm

    bstar = sapply(b,function(x) -rowSums(x*d))

    #Build the cross-derivative block in the big Jacobian A
    a = matrix(NA,ncol=m,nrow=s)
    a[1,] = colSums(bstar)
    a[2:(t+1),] = t(apply(as.matrix(xhat), 2, function(x) colSums(x * bstar)))
    a[(t+2):(t+q+1),] = t(apply(as.matrix(W.main.std), 2, function(x) colSums(x * bstar))) # A2


    w = diag(p)

    #∑pi(pi)[1,xhat,Wi][1,xhat,Wi]^T
    Astar  =  t(cbind(1,xhat, W.main.std)) %*% w %*% cbind(1,xhat, W.main.std) /nm #A3

    # assemble the entire big matrix A
    A =  rbind(cbind(diag(rep(-1,m)),matrix(0,nrow=m,ncol=s)),#A1
               cbind(a,Astar))
    A = rbind(cbind(diag(rep(-1,t+q),nrow=t+q),matrix(0,nrow=t+q,ncol=m+s)),
              cbind(matrix(0,nrow=m+s,ncol=t+q),A))

    ##mus
    mu  = c(muz_std, muw_std)
    dmgi = sapply(1:(t+q),function(x) (Z[x,]-mu[x])/nm) # first line of g(u)
    ###sigmas
    ###diag
    dgi = t(((Z-mu)^2-diag(sigmaz))/nm)    ###off-diag
    ogi = sapply(1:(t+q-1),function(x) sapply((x+1):(t+q), function(y) # thrid line of g(u)
      ((Z[x,]-mu[x])*(Z[y,]-mu[y])-sigmaz[x,y])/nm))

    ###within variance
    ###diag
    zbar = sapply(z.rep.std,function(y) rowMeans(y,na.rm = T))
    r0 = r[indicator==0]

    dgi.0 = if(t==1){
      sapply(1:t,function(y) (rowSums((z.rep.std[[y]]-zbar[,y])^2,na.rm = T)-sigma*(r0-1))/sum(r0-1))
    }else sapply(1:t,function(y) (rowSums((z.rep.std[[y]]-zbar[,y])^2,na.rm = T)-diag(sigma)[y]*(r0-1))/sum(r0-1))
    ###off-diag

    if(t>1){ogi.0 = sapply(1:(t-1),function(x) sapply((x+1):t, function(y)
      (rowSums((z.rep.std[[x]]-zbar[,x])*(z.rep.std[[y]]-zbar[,y]),na.rm = T)-sigma[x,y]*(r0-1))/sum(r0-1)))
    }

    ###betas
    # last three line of g(u)
    gi.beta = cbind(1,xhat,W.main.std)*(Y-p)/nm


    gi.1 = cbind(dmgi,dgi,matrix(unlist(ogi),nrow=nm),gi.beta)
    B1 = t(gi.1)%*%gi.1


    gi.2 = if(t==1){
      dgi.0
    } else cbind(dgi.0,matrix(unlist(ogi.0),nrow=nr))

    B2 = t(gi.2)%*%gi.2

    m1 = (t+q)*(t+q+1)/2
    m2 = t*(t+1)/2

    #Building full B matrix
    B = rbind(cbind(B1[1:(t+q+m1),1:(t+q+m1)],matrix(0,t+q+m1,m2),B1[1:(t+q+m1),(t+q+m1+1):(t+q+m1+s)]),
              cbind(matrix(0,m2,t+q+m1),B2,matrix(0,m2,s)),
              cbind(B1[(t+q+m1+1):(t+q+m1+s),1:(t+q+m1)],matrix(0,s,m2),B1[(t+q+m1+1):(t+q+m1+s),(t+q+m1+1):(t+q+m1+s)]))

    # Sandwich formula
    cov2 = solve(A)%*%B%*%t(solve(A))
    var3 = cov2[(t+q+m+1):(t+q+m+s),(t+q+m+1):(t+q+m+s)]

    tab3 = summary(fit2)$coefficients
    tab3[,2] = sqrt(diag(cov2)[(t+q+m+1):(t+q+m+s)])
    tab3[,1:2] = tab3[,1:2]/c(1,sdz,sdw)
    CI.low = tab3[,1]-1.96*tab3[,2]
    CI.high = tab3[,1]+1.96*tab3[,2]
    tab3[,3] = tab3[,1]/tab3[,2]
    tab3[,4] = 2*pnorm(tab3[,3],lower.tail=FALSE)
    tab3 = cbind(tab3,exp(cbind(OR = tab3[, 1],CI.low,CI.high)))


    rownames(tab3) = sub("^xhat", "", rownames(tab3))
    rownames(tab3) = sub("^W.main.std", "", rownames(tab3))

    return(list(
      `Sandwich Corrected estimates` = tab3
    ))

  }

}
