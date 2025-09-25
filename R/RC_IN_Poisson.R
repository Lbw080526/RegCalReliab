#' Regression Calibration (Internal) for Poisson Models via a Formula Interface
#'
#' @description
#' A convenient wrapper that parses a measurement-error formula, prepares data,
#' runs the naive Poisson log-linear model, performs regression calibration, and
#' applies the sandwich correction that accounts for estimating the measurement model,
#' using internal replicate data.
#'
#' @param formula A formula or character string like
#'   "Y ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + weight".
#'   Terms of the form var(rep1, rep2, ...) are treated as error-prone exposures
#'   with replicates found in main_data. Ordinary terms are treated as covariates W.
#' @param main_data Data frame holding the outcome, replicate error-prone exposures,
#'   and any non-error covariates.
#' @param link Character string specifying the link function for the Poisson model
#'   (default = "log").
#' @param return_details Logical; if TRUE, returns additional internals (xhat, icc, etc.).
#'
#' @return A list with two tidy tables:
#' \itemize{
#'   \item \code{uncorrected}: naive Poisson regression estimates.
#'   \item \code{corrected}: sandwich-corrected RC estimates (your final results).
#' }
#' If \code{return_details = TRUE}, also returns intermediate objects.
#'
#' @examples
#' RC_IN_Poisson(
#'   formula = Y ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + weight,
#'   main_data = main,
#'   link = "log"
#' )
#'
#' @noRd
#' @export

RC_IN_Poisson <- function(formula,
                          main_data,
                          link = "log",
                          return_details = FALSE) {

  # ---- 0) Validate family ----
  if (!is.character(link) || length(link) != 1) {
    stop("`link` must be a single character string such as 'log', 'identity', or 'sqrt'.")
  }
  # Build the Poisson family from the link
  family = poisson(link = link)

  # ---- 1) Normalize formula to string and parse ----
  formula_str <- if (inherits(formula, "formula")) {
    paste(deparse(formula), collapse = "")
  } else {
    as.character(formula)
  }

  parsed <- parse_and_extract_internal(
    main_data   = main_data,
    formula_str = formula_str
  )
  # parsed returns: r, z, W, Y

  # ---- 2) Prepare & standardize ----
  prep <- prepare_data_in(
    r = parsed$r,
    z = parsed$z,
    W = parsed$W,
    Y = parsed$Y
  )
  # prep returns: zbar, z.std, W.std, sds, means, Y, r

  # ---- 3) Naive Poisson regression ----
  naive <- naive_analysis_in_poisson(
    Y     = prep$Y,
    zbar  = prep$zbar,
    W.std = prep$W.std,
    sdz   = prep$sds[["z"]],
    sdw   = prep$sds[["w"]]
  )

  # ---- 4) Regression calibration ----
  rc <- reg_calibration_in_poisson(
    Y     = prep$Y,
    var1  = naive$var1,
    zbar  = prep$zbar,
    z.std = prep$z.std,
    W.std = prep$W.std,
    muz   = prep$means[["z"]],
    muw   = prep$means[["w"]],
    sdz   = prep$sds[["z"]],
    sdw   = prep$sds[["w"]],
    r     = prep$r
  )

  # ---- 5) Sandwich variance estimation ----
  sand <- sandwich_estimator_in_poisson(
    xhat        = rc$xhat,
    zbar        = prep$zbar,
    z.std       = prep$z.std,
    r           = prep$r,
    Y           = prep$Y,
    v12star     = rc$v12star,
    beta.fit2   = rc$beta.fit2,
    W.std       = prep$W.std,
    sigma       = rc$sigma,
    sigmawithin = rc$sigmawithin,
    sigmazstar  = rc$sigmazstar,
    sigmazhat   = rc$sigmazhat,
    sdz         = prep$sds[["z"]],
    sdw         = prep$sds[["w"]],
    muz         = prep$means[["z"]],
    muw         = prep$means[["w"]],
    fit2        = rc$fit2,
    v           = rc$v
  )

  out <- list(
    uncorrected = naive[["Naive estimates"]],
    corrected   = sand[["Sandwich Corrected estimates"]]
  )

  if (return_details) {
    out$details <- list(
      parsed   = parsed,
      prepared = prep,
      rc       = rc
    )
  }

  class(out) <- c("RC_IN_poisson_result", class(out))
  out
}
