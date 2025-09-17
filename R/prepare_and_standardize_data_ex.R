#' Prepare and Standardize Data (External Reliability Study)
#'
#' This function prepares and standardizes the data for the regression calibration process.
#'
#' @param z.main Matrix (n_m x t) of error-prone exposures for the main study.
#' @param r Integer vector of replicate counts for reliability-study subjects.
#' @param z.rep List of length t; each element is an (n_r x r_i) matrix of repeated measurements
#'   for the corresponding exposure in the reliability study.
#' @param W Optional matrix of covariates for the main study.
#' @param Y Binary outcome vector for the main-study subjects.
#'
#' @return A list containing standardized data for the main study and reliability study.
#' @noRd
#' @export


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
