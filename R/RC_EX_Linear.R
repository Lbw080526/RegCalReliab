#' Regression Calibration (External) for Linear Models via a Formula Interface
#'
#' @description
#' Wrapper that parses a measurement-error formula, prepares data,
#' runs the naive linear model, performs regression calibration, and applies
#' the sandwich correction that accounts for estimating the measurement model.
#'
#' @param formula A formula or character string like
#'   "dbp ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + weight".
#'   Terms of the form `var(rep1, rep2, ...)` are treated as error-prone exposures
#'   with replicates found in `rep_data`. Ordinary terms are treated as covariates W.
#' @param main_data Data frame holding the outcome, the main-study error-prone
#'   exposures (the prefixes), and any non-error covariates.
#' @param rep_data  Data frame holding the replicate columns referenced in `formula`.
#' @param family    A gaussian GLM family (must be `gaussian(link = "identity")`
#'   or `gaussian(link = "linear")`).
#' @param return_details Logical; if TRUE, returns additional internals (xhat, icc, etc.).
#'
#' @return A list with two tidy tables:
#' \itemize{
#'   \item \code{uncorrected}: naive linear regression estimates.
#'   \item \code{corrected}: sandwich-corrected RC estimates.
#' }
#' If \code{return_details = TRUE}, also returns intermediate objects.
#'
#' @examples
#' RC_EX_Linear(
#'   formula = dbp ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + weight,
#'   main_data = main, rep_data = ers, family = gaussian("identity")
#' )
#' @noRd

RC_EX_Linear <- function(formula,
                         main_data,
                         rep_data,
                         link = "linear",
                         return_details = FALSE) {

  # ---- 0) Validate link ----
  if (!is.character(link) || length(link) != 1) {
    stop("`link` must be a single character string such as 'linear' or 'identity'.")
  }

  # Treat 'linear' as an alias for 'identity'
  if (link == "linear") link <- "identity"

  # Build the gaussian family with the chosen link
  family <- gaussian(link = link)

  # ---- 1) Normalize formula to string and parse ----
  formula_str <- if (inherits(formula, "formula")) {
    paste(deparse(formula), collapse = "")
  } else {
    paste(as.character(formula), collapse = "")
  }

  parsed <- parse_and_extract(
    main_data     = main_data,
    external_data = rep_data,
    formula_str   = formula_str
  )

  # ---- 2) Prepare & standardize ----
  prep <- prepare_data_ex(
    z.main = parsed$z.main,
    r      = parsed$r,
    z.rep  = parsed$z.rep,
    W      = parsed$W,
    Y      = parsed$Y
  )

  # ---- 3) Naive linear regression ----
  naive <- naive_analysis_ex_linear(
    z.main.std  = prep$z.main.std,
    W.main.std  = prep$W.main.std,
    Y           = prep$Y,
    sdz         = prep$sds[["z"]],
    sdw         = prep$sds[["w"]]
  )

  # ---- 4) Regression calibration ----
  rc <- reg_calibration_ex_linear(
    z.main.std   = prep$z.main.std,
    z.rep.std    = prep$z.rep.std,
    r            = prep$r,
    W.main.std   = prep$W.main.std,
    Y            = prep$Y,
    muz          = prep$means[["z"]],
    muw          = prep$means[["w"]],
    sdz          = prep$sds[["z"]],
    sdw          = prep$sds[["w"]],
    indicator    = prep$indicator
  )

  # ---- 5) Sandwich variance estimation ----
  sand <- sandwich_estimator_ex_linear(
    xhat        = rc$xhat,
    z.main.std  = prep$z.main.std,
    z.rep.std   = prep$z.rep.std,
    r           = prep$r,
    Y           = prep$Y,
    indicator   = prep$indicator,
    v12star     = rc$v12star,
    beta.fit2   = rc$beta.fit2,
    W.main.std  = prep$W.main.std,
    sigma       = rc$sigma,
    sigmaz      = rc$sigmaz,
    sigmawithin = rc$sigmawithin,
    sdz         = prep$sds[["z"]],
    sdw         = prep$sds[["w"]],
    fit2        = rc$fit2
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

  class(out) <- c("RC_EX_linear_result", class(out))
  out
}
