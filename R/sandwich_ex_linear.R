#' Sandwich Variance Estimator for External Logistic Regression Calibration
#'
#' \code{sandwich_estimator_ex_linear()} computes robust (sandwich) standard errors
#' for regression calibration in linear regression using external reliability data.
#' It propagates uncertainty from estimating the calibration (variance components)
#' into the final coefficient variances, yielding valid confidence intervals and
#' p-values even when measurement error is present.
#'
#' @param xhat Numeric matrix of calibrated exposures (\eqn{n_m \times t}),
#'   typically obtained from \code{\link{reg_calibration_ex_linear}}.
#' @param z.main.std Numeric matrix of standardized main-study exposures
#'   (\eqn{n_m \times t}).
#' @param z.rep.std Named list of standardized replicate measurements from the
#'   external reliability study. Each list element is a matrix of dimension
#'   \eqn{n_r \times r_i} for a particular exposure.
#' @param r Integer vector of replicate counts for all subjects
#'   (length \eqn{n_m + n_r}); main-study subjects should have value 1.
#' @param Y Numeric outcome vector of length \eqn{n_m}.
#' @param indicator Binary vector of length \eqn{n_m + n_r} indicating main-study
#'   (1) vs. reliability-study (0) subjects.
#' @param v12star Calibration matrix (\eqn{\Sigma_X \Sigma_Z^{-1}\Sigma_{XZ}})
#'   used to form \code{xhat}.
#' @param beta.fit2 Numeric vector of regression coefficients from the corrected
#'   linear model (\code{fit2}).
#' @param W.main.std Optional numeric matrix of standardized error-free
#'   covariates (\eqn{n_m \times q}). If not provided, the model is fit with
#'   exposures only.
#' @param sigma Within-person variance-covariance matrix of the exposures.
#' @param sigmaz Total variance-covariance matrix of the observed exposures
#'   (and covariates, if applicable).
#' @param sigmawithin Estimated within-person covariance matrix (same dimension
#'   as \code{sigmaz}).
#' @param sdz Numeric vector of standard deviations of the unstandardized exposures.
#' @param sdw Optional numeric vector of standard deviations of the unstandardized covariates.
#' @param fit2 The \code{lm} object returned by
#'   \code{\link{reg_calibration_ex_linear}}.
#'
#' @return A list with one component:
#' \describe{
#'   \item{\code{Sandwich Corrected estimates}}{Matrix of regression calibration
#'         estimates with sandwich-corrected standard errors, t-values,
#'         p-values, and 95\% confidence intervals on the original scale.}
#' }
#'
#' @details
#' This function implements the three-block estimating equation described in
#' Carroll et al. (2006) and White (1980):
#' \enumerate{
#'   \item \strong{Variance components}: estimates total and within-person covariance.
#'   \item \strong{Calibration}: computes \eqn{X^\text{hat}} using the
#'         estimated variance components.
#'   \item \strong{Outcome model}: fits a linear regression and applies a
#'         sandwich (robust) variance formula to combine the uncertainty from
#'         all stages.
#' }
#'
#' @examples
#' set.seed(1)
#' # Simulated main-study data: 80 subjects, 1 exposure
#' z <- matrix(rnorm(80), ncol = 1)
#' colnames(z) <- "sbp"
#' Y <- 2 + 0.5 * z + rnorm(80)
#'
#' # Reliability study: 40 subjects, 2 replicates
#' z.rep <- list(sbp = matrix(rnorm(40 * 2), nrow = 40))
#' zbar  <- sapply(z.rep, rowMeans)
#'
#' # Standardize data
#' sdz <- apply(z, 2, sd)
#' z.main.std <- scale(z)
#' z.rep.std  <- list(sbp = scale(z.rep$sbp))
#' r <- c(rep(1, 80), rep(2, 40))
#' indicator <- c(rep(1, 80), rep(0, 40))
#'
#' # (In practice, run reg_calibration_ex_linear() first to get xhat, fit2, etc.)
#' # sandwich_estimator_ex_linear(xhat, z.main.std, z.rep.std, r, Y,
#' #   indicator, v12star, beta.fit2, sigma, sigmaz, sigmawithin,
#' #   sdz, fit2 = fit2)
#'
#' @references
#' Carroll RJ, Ruppert D, Stefanski LA, Crainiceanu CM. \emph{Measurement Error
#' in Nonlinear Models: A Modern Perspective}. Chapman & Hall/CRC, 2006.
#' White H. A heteroskedasticity-consistent covariance matrix estimator and a
#' direct test for heteroskedasticity. \emph{Econometrica}. 1980;48(4):817–838.
#'
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

  list(`Sandwich Corrected estimates` = tab3)
}
