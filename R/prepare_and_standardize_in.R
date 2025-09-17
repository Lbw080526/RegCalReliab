#' Prepare and Standardize Data (Internal Replicate Study)
#'
#' Combines replicate measurements \code{z} from an internal reliability study,
#' constructs replicate counts \code{r}, computes subject-level averages
#' \code{zbar}, and standardizes exposures (and optional covariates \code{W})
#' using means/SDs across all replicates. Returns standardized replicate-level
#' and subject-level data for subsequent naive, calibration, and sandwich analyses.
#'
#' @param r Integer vector of replicate counts (length n).
#' @param z List of length \code{t}; each element is an (n x r_i) matrix of
#'   repeated measurements for one error-prone exposure. Elements may differ
#'   in number of columns (replicates).
#' @param W Optional matrix of additional covariates; default \code{NULL}.
#' @param Y Outcome vector (length n).
#'
#' @return A list with:
#' \item{z.std}{List of length \code{t}; each element is an (n x r_max) matrix
#'   of standardized replicates, padded with \code{NA} if fewer than the maximum
#'   number of replicates.}
#' \item{zbar}{Matrix (n x t) of standardized subject-level averages of replicates.}
#' \item{W.std}{If \code{W} supplied: standardized covariate matrix (n x q).}
#' \item{Y}{Outcome vector (length n).}
#' \item{r}{Integer vector of replicate counts (length n).}
#' \item{means}{List of means used for standardization: \code{z}, and \code{w}
#'   if covariates are provided.}
#' \item{sds}{List of standard deviations used for standardization: \code{z}, and
#'   \code{w} if covariates are provided.}
#'
#' @details
#' Standardization is performed using the pooled mean and
#' standard deviation across all replicates. Replicate matrices \code{z}
#' are padded with \code{NA} to achieve equal width
#' across exposures. Subject-level averages \code{zbar}
#' are computed as row means of standardized replicates, ignoring \code{NA}.
#'
#' @noRd
#' @export



prepare_data_in = function(r, z, W = NULL, Y) {


  n = length(r)
  t = length(z)

  z_name = names(z)

  z = lapply(1:t, function(x) {
    zx <- z[[x]]
    pad_cols <- max(r) - ncol(zx)
    if (pad_cols > 0) {
      zx <- cbind(zx, matrix(NA, nrow = n, ncol = pad_cols))
    }
    zx
  })


  meanz = sapply(z, function(x) colMeans(x, na.rm = TRUE))[1, ]
  sdz = sapply(z, function(x) apply(x, 2, function(y) sd(y, na.rm = TRUE)))[1, ]
  z.std = sapply(1:t, function(x) (z[[x]] - meanz[x]) / sdz[x], simplify = FALSE)
  names(z.std) = z_name

  if (is.null(W)) {
    zbar = sapply(z.std, function(y) rowMeans(y, na.rm = TRUE))
    names(zbar) = z_name


    return(list(
      means = list(z = meanz),
      z.std = z.std,
      sds   = list(z = sdz),
      Y     = Y,
      zbar  = zbar,
      r = r
    ))
  }
  else {
    W = as.matrix(W)
    t = length(z)
    q = ncol(W)
    meanw = colMeans(W)
    sdw = apply(W, 2, function(y) sd(y, na.rm = TRUE))
    W.new = sapply(1:q, function(x) (W[, x] - meanw[x]) / sdw[x])
    colnames(W.new) = colnames(W)
    zbar = sapply(z.std, function(y) rowMeans(y, na.rm = TRUE))
    names(zbar) = names(z)

    return(list(
      W.std = W.new,
      means = list(z = meanz, w = meanw),
      sds   = list(z = sdz, w = sdw),
      zbar  = zbar,
      Y = Y,
      z.std = z.std,
      r = r
    ))
  }
}
