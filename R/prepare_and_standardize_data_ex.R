#' Prepare and Standardize Data for External Reliability Study
#'
#' \code{prepare_data_ex()} merges main-study and external-reliability data
#' into a unified structure for regression calibration and standardizes all
#' error-prone exposures (and optionally covariates).
#' This function takes the main-study error-prone measurements,
#' their repeated measurements from the external reliability study,
#' the number of replicates, and optional error-free covariates,
#' and returns standardized matrices and associated information.
#'
#' Typical usage:
#' \enumerate{
#'   \item Pads each main-study exposure to match the maximum number of replicates.
#'   \item Concatenates main-study and reliability-study observations into a single dataset.
#'   \item Computes means and standard deviations for each exposure (and covariates if supplied).
#'   \item Returns standardized main-study exposures, replicate measurements,
#'         replicate counts, and other auxiliary information.
#' }
#'
#' @param z.main A numeric matrix of dimension \eqn{n_m \times t}, where
#'   \eqn{n_m} is the main-study sample size and \eqn{t} is the number of
#'   error-prone exposures. Each column is a single main-study measurement
#'   for that exposure.
#' @param r Integer vector of length \eqn{n_r} giving the number of replicates
#'   available for each subject in the external reliability study.
#' @param z.rep A named list of length \eqn{t}; each element is a matrix of
#'   dimension \eqn{n_r \times r_i} containing replicate measurements of the
#'   corresponding exposure in the reliability study. The list names must match
#'   the column names of \code{z.main}.
#' @param W Optional numeric matrix of error-free covariates for the main study
#'   (\eqn{n_m \times q}). If supplied, all covariates are centered and scaled.
#' @param Y Numeric (0/1) vector of length \eqn{n_m} giving the binary outcome
#'   for the main study.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{z.main.std}}{Standardized main-study exposures (\eqn{n_m \times t}).}
#'   \item{\code{z.rep.std}}{List of standardized replicate measurements, each a matrix of dimension \eqn{n_r \times r_i}.}
#'   \item{\code{r}}{Integer vector of replicate counts of length \eqn{n_m + n_r}
#'         (main-study subjects are coded as 1).}
#'   \item{\code{W.main.std}}{Standardized main-study covariates (\eqn{n_m \times q}) if \code{W} is supplied.}
#'   \item{\code{Y}}{Binary outcome vector of length \eqn{n_m}.}
#'   \item{\code{indicator}}{Binary vector of length \eqn{n_m + n_r} identifying
#'         main-study subjects (1) vs. reliability-study subjects (0).}
#'   \item{\code{means}}{List of means used for standardization: \code{z} for exposures,
#'         and \code{w} for covariates if applicable.}
#'   \item{\code{sds}}{List of standard deviations used for standardization: \code{z} for exposures,
#'         and \code{w} for covariates if applicable.}
#' }
#'
#' @examples
#' set.seed(123)
#' # Main study: 6 subjects, 2 error-prone exposures
#' z.main <- matrix(rnorm(6 * 2), nrow = 6, ncol = 2)
#' colnames(z.main) <- c("sbp", "chol")
#'
#' # Reliability study: 4 subjects, each with 2 replicates
#' z.rep <- list(
#'   sbp  = matrix(rnorm(4 * 2), nrow = 4),
#'   chol = matrix(rnorm(4 * 2), nrow = 4)
#' )
#'
#' # Replicate counts
#' r <- rep(2, 4)
#'
#' # Optional covariates and binary outcome
#' W <- matrix(rnorm(6 * 2), nrow = 6)
#' colnames(W) <- c("age", "sex")
#' Y <- rbinom(6, 1, 0.5)
#'
#' # Prepare and standardize data
#' prep <- prepare_data_ex(z.main = z.main,
#'                         r = r,
#'                         z.rep = z.rep,
#'                         W = W,
#'                         Y = Y)
#' str(prep)
#'
#' @noRd


prepare_data_ex <- function(z.main, r, z.rep, W = NULL, Y) {
  nm = nrow(z.main)             # main-study sample size
  nr = length(r)                # reliability-study sample size
  r  = c(rep(1, nm), r)         # unify replicate counts: main=1, reliability = r
  n = nm + nr

  var_names = names(z.rep)
  t = length(var_names)

  max_reps = max(sapply(z.rep, ncol))

  # Build combined data
  z = lapply(var_names, function(vname) {
    main_padded = cbind(z.main[, vname], matrix(NA, nrow = nm, ncol = max_reps - 1))
    reliability_data = z.rep[[vname]]
    if (ncol(reliability_data) < max_reps) {
      reliability_data = cbind(reliability_data, matrix(NA, nrow = nr, ncol = max_reps - ncol(reliability_data)))
    }
    rbind(main_padded, reliability_data)
  })
  names(z) = var_names

  indicator = c(rep(1, nm), rep(0, nr))

  # Means and SDs
  meanz = sapply(z, function(x) colMeans(x, na.rm = TRUE))[1, ]
  sdz   = sapply(z, function(x) apply(x, 2, sd, na.rm = TRUE))[1, ]

  # Standardize z
  z.new = mapply(function(mat, mu, sdv) (mat - mu) / sdv,
                 z, meanz, sdz, SIMPLIFY = FALSE)

  z.main.std = sapply(z.new, function(y) y[indicator == 1, 1, drop = FALSE])
  colnames(z.main.std) = var_names
  z.rep.std = lapply(z.new, function(y) y[indicator == 0, ])

  if (is.null(W)) {
    return(list(
      z.main.std = z.main.std,
      z.rep.std  = z.rep.std,
      r          = r,
      Y          = Y,
      indicator  = indicator,
      means      = list(z = meanz),
      sds        = list(z = sdz)
    ))
  }

  # If W exists
  q = ncol(W)
  W = as.matrix(W)
  meanw = colMeans(W)
  sdw   = apply(W, 2, sd)
  W.new = scale(W, center = meanw, scale = sdw)
  W.main.std = matrix(W.new, nrow = nm)
  colnames(W.main.std) = colnames(W)

  return(list(
    z.main.std = z.main.std,
    z.rep.std  = z.rep.std,
    r          = r,
    W.main.std = W.main.std,
    Y          = Y,
    indicator  = indicator,
    means      = list(z = meanz, w = meanw),
    sds        = list(z = sdz, w = sdw)
  ))
}
