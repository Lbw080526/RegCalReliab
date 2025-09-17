#' Sandwich Variance Estimator for External Linear Regression Calibration
#'
#' Robust (sandwich) SEs for regression-calibrated **linear** models with an
#' external reliability study. Works with or without additional covariates W.
#'
#' @param xhat        Matrix of calibrated predictors for the main study (n_m x t).
#' @param z.main.std  Standardized main-study exposures (n_m x t).
#' @param z.rep.std   List of standardized replicate matrices for reliability data.
#' @param r           Integer vector of replicate counts (length n_r).
#' @param Y           Numeric outcome vector (length n_m).
#' @param indicator   Binary vector: 1 = main-study rows, 0 = reliability rows (length n_m + n_r).
#' @param v12star     Between-person exposure block used to form xhat.
#' @param beta.fit2   Coefficients from the corrected **linear** model (fit on xhat [+ W]).
#' @param W.main.std  Optional standardized covariates (n_m x q).
#' @param sigma       Within-person covariance (t x t).
#' @param sigmaz      Total covariance of (Z[,W]) used in calibration: (t x t) or (t+q) x (t+q).
#' @param sigmawithin Within-person covariance (same block structure as sigmaz).
#' @param sdz         Vector of SDs used to standardize exposures (length t).
#' @param sdw         Optional vector of SDs used to standardize covariates (length q).
#' @param fit2        The fitted linear model object from the calibration step (for coefficient table).
#'
#' @return list with element `Sandwich Corrected estimates`: coefficient table with
#'         Estimate, Std. Error, z value, Pr(>|z|), CI.low, CI.high,
#'         transformed back to original exposure/covariate scales via (sdz, sdw).
#' @noRd
#' @export

sandwich_estimator_ex_linear <- function(
    xhat, z.main.std, z.rep.std, r, Y, indicator, v12star, beta.fit2,
    W.main.std = NULL, sigma, sigmaz, sigmawithin, sdz, sdw = NULL, fit2
) {
  ## ---------- 0) Dimensions & helpers ----------
  nm <- nrow(z.main.std)
  t  <- ncol(z.main.std)
  xhat <- as.matrix(if (is.null(dim(xhat))) matrix(xhat, ncol = t) else xhat)

  # Means used in g(u): MUST be on the same scale as Z below (standardized)
  muz_std <- as.numeric(colMeans(z.main.std))
  if (!is.null(W.main.std)) muw_std <- as.numeric(colMeans(W.main.std))

  sigmaz.inv <- solve(sigmaz)   # linear algebra matches reference implementation

  ## ======================================================================
  ## CASE 1: NO W
  ## ======================================================================
  if (is.null(W.main.std)) {
    # Model mean & residuals under linear model
    X_lin <- cbind(1, xhat)              # design used in corrected linear fit
    mu    <- as.vector(X_lin %*% beta.fit2)
    resid <- Y - mu

    # Z matrix for covariance derivatives (transpose so rows = variables)
    Z <- t(z.main.std)                   # (t x nm)

    # --- derivative blocks wrt covariance params ---
    m0 <- matrix(0, t, t)
    c0 <- matrix(0, t, t)

    ## diag blocks for Sigma_Z (top-left)
    ddv.x <- lapply(1:t, function(x){ a <- rep(0,t); a[x] <- 1; diag(a) })
    ddm.x <- lapply(1:t, function(x){ a <- rep(0,t); a[x] <- 1; diag(a) })

    db.x <- lapply(1:t, function(x) {
      t((ddv.x[[x]] %*% sigmaz.inv -
           v12star %*% sigmaz.inv %*% ddm.x[[x]] %*% sigmaz.inv) %*% Z)
    })

    ## off-diag for Sigma_Z
    if (t > 1) {
      odv.x <- lapply(1:(t-1), function(x) lapply((x+1):t, function(y){
        c <- c0; c[x,y] <- c[y,x] <- 1; c
      }))
      odm.x <- lapply(1:(t-1), function(x) lapply((x+1):t, function(y){
        m <- m0; m[x,y] <- m[y,x] <- 1; m
      }))
      ob.x <- lapply(1:(t-1), function(x) lapply(1:(t-x), function(y)
        t((odv.x[[x]][[y]] %*% sigmaz.inv -
             v12star %*% sigmaz.inv %*% odm.x[[x]][[y]] %*% sigmaz.inv) %*% Z)))
    }

    ## within-person sigma (diagonal block)
    db.0 <- lapply(1:t, function(x) t((ddv.x[[x]] %*% sigmaz.inv) %*% Z))
    if (t > 1) {
      ob.0 <- lapply(1:(t-1), function(x) lapply(1:(t-x), function(y)
        t((odv.x[[x]][[y]] %*% sigmaz.inv) %*% Z)))
    }

    m <- t*(t+1)/2 + t*(t+1)/2          # number of covariance params
    s <- t + 1                          # number of beta params

    # Pack derivative blocks into list b of nm x t matrices
    blocks_list <- if (t > 1) list(db.x, ob.x, db.0, ob.0) else list(db.x, db.0)
    B_un <- matrix(unlist(blocks_list), nrow = nm)
    b <- lapply(seq_len(m), function(ix) {
      M <- B_un[, (t*ix - t + 1):(t*ix), drop = FALSE]
      if (is.null(dim(M))) matrix(M, nrow = nm, ncol = t) else M
    })

    # Linear case: d = beta_Z / nm (no p(1-p) factor here)
    d_linear <- matrix(rep(beta.fit2[2:(t+1)], each = nm), nrow = nm, ncol = t) / nm
    bstar    <- sapply(b, function(mat) -rowSums(mat * d_linear))

    # Cross-derivative block A12 (a)
    a <- matrix(NA, nrow = s, ncol = m)
    a[1, ] <- colSums(bstar)
    a[2:(t+1), ] <- t(apply(xhat, 2, function(xx) colSums(xx * bstar)))

    # Sensitivity wrt beta (A22)
    Astar <- crossprod(X_lin) / nm

    # Build big A: rows/cols are [mu_z params | cov params | beta]
    A <- rbind(cbind(diag(rep(-1, m)), matrix(0, m, s)),
               cbind(a, Astar))
    A <- rbind(cbind(diag(rep(-1, t), t), matrix(0, t, m + s)),
               cbind(matrix(0, m + s, t), A))

    ## ---- Score pieces g(u) for B ----
    # main-study part
    dmgi <- sapply(1:t, function(x) (Z[x,] - muz_std[x]) / nm)
    dgi  <- t(((Z - muz_std)^2 - diag(sigmaz)) / nm)
    if (t > 1) {
      ogi <- lapply(1:(t-1), function(x) sapply((x+1):t, function(y)
        ((Z[x,]-muz_std[x]) * (Z[y,]-muz_std[y]) - sigmaz[x,y]) / nm))
    }

    gi.beta <- X_lin * resid / nm

    gi.1 <- if (t == 1) cbind(dmgi, dgi, gi.beta) else
      cbind(dmgi, dgi, matrix(unlist(ogi), nrow = nm), gi.beta)
    B1 <- crossprod(gi.1)

    # reliability part
    zbar <- sapply(z.rep.std, function(y) rowMeans(y, na.rm = TRUE))
    r0   <- r[indicator == 0]
    dgi.0 <- if (t == 1) {
      sapply(1:t, function(y)
        (rowSums((z.rep.std[[y]] - zbar[, y])^2, na.rm = TRUE) - sigma * (r0 - 1)) / sum(r0 - 1))
    } else {
      sapply(1:t, function(y)
        (rowSums((z.rep.std[[y]] - zbar[, y])^2, na.rm = TRUE) - diag(sigma)[y] * (r0 - 1)) / sum(r0 - 1))
    }
    if (t > 1) {
      ogi.0 <- lapply(1:(t-1), function(x) sapply((x+1):t, function(y)
        (rowSums((z.rep.std[[x]] - zbar[, x]) * (z.rep.std[[y]] - zbar[, y]), na.rm = TRUE) -
           sigma[x,y] * (r0 - 1)) / sum(r0 - 1)))
    }

    gi.2 <- if (t == 1) dgi.0 else cbind(dgi.0, matrix(unlist(ogi.0), nrow = length(r0)))
    B2   <- crossprod(gi.2)

    m1 <- t*(t+1)/2
    m2 <- t*(t+1)/2

    B <- rbind(
      cbind(B1[1:(t+m1), 1:(t+m1)], matrix(0, t+m1, m2), B1[1:(t+m1), (t+m1+1):(t+m1+s)]),
      cbind(matrix(0, m2, t+m1), B2, matrix(0, m2, s)),
      cbind(B1[(t+m1+1):(t+m1+s), 1:(t+m1)], matrix(0, s, m2),
            B1[(t+m1+1):(t+m1+s), (t+m1+1):(t+m1+s)])
    )

    cov2 <- solve(A) %*% B %*% t(solve(A))
    idx_beta <- (t + m + 1):(t + m + s)

    tab3 <- summary(fit2)$coefficients
    tab3[, 2] <- sqrt(diag(cov2)[idx_beta])
    tab3[, 1:2] <- tab3[, 1:2] / c(1, sdz)
    tab3[, 3] <- tab3[, 1] / tab3[, 2]
    tab3[, 4] <- 2 * pnorm(tab3[, 3], lower.tail = FALSE)
    CI.low  <- tab3[, 1] - 1.96 * tab3[, 2]
    CI.high <- tab3[, 1] + 1.96 * tab3[, 2]
    tab3 <- cbind(tab3, CI.low = CI.low, CI.high = CI.high)
    rownames(tab3) <- sub("^xhat", "", rownames(tab3))

    return(list(`Sandwich Corrected estimates` = tab3))
  }

  ## ======================================================================
  ## CASE 2: WITH W
  ## ======================================================================
  colnames(xhat)       <- colnames(z.main.std)
  colnames(W.main.std) <- colnames(W.main.std)
  nm <- nrow(z.main.std)
  q  <- ncol(W.main.std)

  X_lin <- cbind(1, xhat, W.main.std)
  mu    <- as.vector(X_lin %*% beta.fit2)
  resid <- Y - mu

  Z <- t(cbind(z.main.std, W.main.std))         # (t+q) x nm

  m0 <- matrix(0, t+q, t+q)
  c0 <- matrix(0, t,   t+q)

  ## Sigma_Z diag
  ddv.x <- lapply(1:t, function(x){
    a <- rep(0,t); a[x] <- 1
    cbind(diag(a), matrix(0, nrow = t, ncol = q))
  })
  ddm.x <- lapply(1:t, function(x){
    a <- rep(0, t+q); a[x] <- 1; diag(a)
  })
  db.x <- lapply(1:t, function(x)
    t((ddv.x[[x]] %*% sigmaz.inv -
         v12star %*% sigmaz.inv %*% ddm.x[[x]] %*% sigmaz.inv) %*% Z))

  ## Sigma_Z off-diag
  if (t > 1) {
    odv.x <- lapply(1:(t-1), function(x) lapply((x+1):t, function(y){
      c <- c0; c[x,y] <- c[y,x] <- 1; c
    }))
    odm.x <- lapply(1:(t-1), function(x) lapply((x+1):t, function(y){
      m <- m0; m[x,y] <- m[y,x] <- 1; m
    }))
    ob.x <- lapply(1:(t-1), function(x) lapply(1:(t-x), function(y)
      t((odv.x[[x]][[y]] %*% sigmaz.inv -
           v12star %*% sigmaz.inv %*% odm.x[[x]][[y]] %*% sigmaz.inv) %*% Z)))
  }

  ## Sigma_ZW (all off-diag)
  odv.xw <- lapply(1:t, function(x) lapply((t+1):(t+q), function(y){
    c <- c0; c[x,y] <- 1; c
  }))
  odm.xw <- lapply(1:t, function(x) lapply((t+1):(t+q), function(y){
    m <- m0; m[x,y] <- m[y,x] <- 1; m
  }))
  ob.xw <- lapply(1:t, function(x) lapply(1:q, function(y)
    t((odv.xw[[x]][[y]] %*% sigmaz.inv -
         v12star %*% sigmaz.inv %*% odm.xw[[x]][[y]] %*% sigmaz.inv) %*% Z)))

  ## Sigma_W
  ddm.w <- lapply((t+1):(t+q), function(x){
    a <- rep(0, t+q); a[x] <- 1; diag(a)
  })
  db.w <- lapply(1:q, function(x)
    t((-v12star %*% sigmaz.inv %*% ddm.w[[x]] %*% sigmaz.inv) %*% Z))
  if (q > 1) {
    odm.w <- lapply((t+1):(t+q-1), function(x) lapply((x+1):(t+q), function(y){
      m <- m0; m[x,y] <- m[y,x] <- 1; m
    }))
    ob.w <- lapply(1:(q-1), function(x) lapply(1:(q-x), function(y)
      t((-v12star %*% sigmaz.inv %*% odm.w[[x]][[y]] %*% sigmaz.inv) %*% Z)))
  }

  ## sigma (within)
  db.0 <- lapply(1:t, function(x) t((-ddv.x[[x]] %*% sigmaz.inv) %*% Z))
  if (t > 1) {
    ob.0 <- lapply(1:(t-1), function(x) lapply(1:(t-x), function(y)
      t((-odv.x[[x]][[y]] %*% sigmaz.inv) %*% Z)))
  }

  m <- (t+q)*(t+q+1)/2 + t*(t+1)/2
  s <- t + q + 1

  blocks_list <- if (t > 1 & q > 1) {
    list(db.x, db.w, ob.x, ob.xw, ob.w, db.0, ob.0)
  } else if (q > 1) {
    list(db.x, db.w, ob.xw, ob.w, db.0)
  } else if (t > 1) {
    list(db.x, db.w, ob.x, ob.xw, db.0, ob.0)
  } else {
    list(db.x, db.w, ob.xw, db.0)
  }
  B_un <- matrix(unlist(blocks_list), nrow = nm)
  b <- lapply(seq_len(m), function(ix){
    M <- B_un[, (t*ix - t + 1):(t*ix), drop = FALSE]
    if (is.null(dim(M))) matrix(M, nrow = nm, ncol = t) else M
  })

  d_linear <- matrix(rep(beta.fit2[2:(t+1)], each = nm), nrow = nm, ncol = t) / nm
  bstar    <- sapply(b, function(M) -rowSums(M * d_linear))

  a <- matrix(NA, nrow = s, ncol = m)
  a[1, ] <- colSums(bstar)
  a[2:(t+1), ] <- t(apply(xhat, 2, function(xx) colSums(xx * bstar)))
  a[(t+2):(t+q+1), ] <- t(apply(W.main.std, 2, function(xx) colSums(xx * bstar)))

  # Sensitivity wrt beta
  Astar <- crossprod(X_lin) / nm
  A <- rbind(cbind(diag(rep(-1, m)), matrix(0, m, s)),
             cbind(a, Astar))
  A <- rbind(cbind(diag(rep(-1, t+q), t+q), matrix(0, t+q, m + s)),
             cbind(matrix(0, m + s, t+q), A))

  ## ---- Score pieces B ----
  mu_all <- c(muz_std, muw_std)
  dmgi <- sapply(1:(t+q), function(x) (Z[x,] - mu_all[x]) / nm)
  dgi  <- t(((Z - mu_all)^2 - diag(sigmaz)) / nm)
  ogi  <- lapply(1:(t+q-1), function(x) sapply((x+1):(t+q), function(y)
    ((Z[x,]-mu_all[x]) * (Z[y,]-mu_all[y]) - sigmaz[x,y]) / nm))

  zbar <- sapply(z.rep.std, function(y) rowMeans(y, na.rm = TRUE))
  r0   <- r[indicator == 0]
  dgi.0 <- if (t == 1) {
    sapply(1:t, function(y)
      (rowSums((z.rep.std[[y]] - zbar[, y])^2, na.rm = TRUE) - sigma * (r0 - 1)) / sum(r0 - 1))
  } else {
    sapply(1:t, function(y)
      (rowSums((z.rep.std[[y]] - zbar[, y])^2, na.rm = TRUE) - diag(sigma)[y] * (r0 - 1)) / sum(r0 - 1))
  }
  if (t > 1) {
    ogi.0 <- lapply(1:(t-1), function(x) sapply((x+1):t, function(y)
      (rowSums((z.rep.std[[x]] - zbar[, x]) * (z.rep.std[[y]] - zbar[, y]), na.rm = TRUE) -
         sigma[x,y] * (r0 - 1)) / sum(r0 - 1)))
  }

  gi.beta <- X_lin * resid / nm

  gi.1 <- cbind(dmgi, dgi, matrix(unlist(ogi), nrow = nm), gi.beta)
  B1   <- crossprod(gi.1)

  gi.2 <- if (t == 1) dgi.0 else cbind(dgi.0, matrix(unlist(ogi.0), nrow = length(r0)))
  B2   <- crossprod(gi.2)

  m1 <- (t+q)*(t+q+1)/2
  m2 <- t*(t+1)/2

  B <- rbind(
    cbind(B1[1:(t+q+m1), 1:(t+q+m1)], matrix(0, t+q+m1, m2),
          B1[1:(t+q+m1), (t+q+m1+1):(t+q+m1+s)]),
    cbind(matrix(0, m2, t+q+m1), B2, matrix(0, m2, s)),
    cbind(B1[(t+q+m1+1):(t+q+m1+s), 1:(t+q+m1)], matrix(0, s, m2),
          B1[(t+q+m1+1):(t+q+m1+s), (t+q+m1+1):(t+q+m1+s)])
  )

  cov2 <- solve(A) %*% B %*% t(solve(A))
  idx_beta <- (t + q + m + 1):(t + q + m + s)

  tab3 <- summary(fit2)$coefficients
  tab3[, 2] <- sqrt(diag(cov2)[idx_beta])
  tab3[, 1:2] <- tab3[, 1:2] / c(1, sdz, sdw)
  tab3[, 3] <- tab3[, 1] / tab3[, 2]
  tab3[, 4] <- 2 * pnorm(tab3[, 3], lower.tail = FALSE)
  CI.low  <- tab3[, 1] - 1.96 * tab3[, 2]
  CI.high <- tab3[, 1] + 1.96 * tab3[, 2]
  tab3 <- cbind(tab3, CI.low = CI.low, CI.high = CI.high)
  rownames(tab3) <- sub("^xhat", "", rownames(tab3))
  rownames(tab3) <- sub("^W\\.main\\.std", "", rownames(tab3))

  list(`Sandwich Corrected estimates` = tab3)
}
