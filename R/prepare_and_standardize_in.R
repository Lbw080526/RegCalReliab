#' Prepare and Standardize Data for Internal Reliability Study
#'
#' \code{prepare_data_in()} processes replicate measurements from an internal
#' reliability study, constructs replicate counts, computes subject-level
#' averages, and standardizes exposures (and optional covariates). It returns
#' replicate-level and subject-level standardized data for downstream naive,
#' calibration, and sandwich analyses.
#'
#' Typical usage:
#' \enumerate{
#'   \item Pads each replicate matrix to match the maximum number of replicates.
#'   \item Computes pooled means and standard deviations across all replicates.
#'   \item Standardizes replicate-level data and subject-level averages.
#'   \item Standardizes optional covariates if supplied.
#'   \item Returns standardized replicates, averages, replicate counts, and
#'         auxiliary information.
#' }
#'
#' @param r Integer vector of length \eqn{n}, giving the number of replicates
#'   available for each subject.
#' @param z A named list of length \eqn{t}; each element is an \eqn{n \times r_i}
#'   matrix of replicate measurements for one error-prone exposure. Different
#'   exposures may have different numbers of replicates.
#' @param W Optional numeric matrix of covariates without measurement error
#'   (\eqn{n \times q}). If supplied, all covariates are centered and scaled.
#' @param Y Outcome vector of length \eqn{n}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{z.std}}{List of standardized replicate matrices, each padded to
#'         the maximum number of replicates.}
#'   \item{\code{zbar}}{Matrix (\eqn{n \times t}) of standardized subject-level
#'         averages of replicates.}
#'   \item{\code{W.std}}{Standardized covariate matrix (\eqn{n \times q}) if
#'         \code{W} is supplied.}
#'   \item{\code{Y}}{Outcome vector of length \eqn{n}.}
#'   \item{\code{r}}{Integer vector of replicate counts of length \eqn{n}.}
#'   \item{\code{means}}{List of means used for standardization: \code{z} for
#'         exposures, and \code{w} for covariates if applicable.}
#'   \item{\code{sds}}{List of standard deviations used for standardization:
#'         \code{z} for exposures, and \code{w} for covariates if applicable.}
#' }
#'
#' @examples
#' set.seed(123)
#' # Internal study: 5 subjects, 1 exposure with 2 replicates
#' z = list(
#'   sbp = cbind(rnorm(5, 120, 15), rnorm(5, 120, 15))
#' )
#'
#' # Replicate counts (each subject has 2 replicates)
#' r = rep(2, 5)
#'
#' # Outcome and optional covariate
#' Y = rbinom(5, 1, 0.5)
#' W = matrix(rnorm(5), ncol = 1)
#' colnames(W) = "age"
#'
#' # Prepare and standardize data
#' prep = prepare_data_in(r = r, z = z, W = W, Y = Y)
#' str(prep)
#'
#' @noRd




prepare_data_in = function(r, z, W = NULL, Y) {


  n = length(r)
  t = length(z)

  z_name = names(z)

  z = lapply(1:t, function(x) {
    zx = z[[x]]
    pad_cols = max(r) - ncol(zx)
    if (pad_cols > 0) {
      zx = cbind(zx, matrix(NA, nrow = n, ncol = pad_cols))
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
