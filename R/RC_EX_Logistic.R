#' Regression Calibration (External) for Logistic Models via a Formula Interface
#'
#' @description
#' A convenient wrapper that parses a measurement-error formula, prepares data,
#' runs the naive logistic model, performs regression calibration, and applies
#' the sandwich correction that accounts for estimating the measurement model.
#'
#' @param formula A formula or character string like
#'   "comb ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + weight".
#'   Terms of the form `var(rep1, rep2, ...)` are treated as error-prone exposures
#'   with replicates found in `rep_data`. Ordinary terms are treated as covariates W.
#' @param main_data Data frame holding the outcome, the main-study error-prone
#'   exposures (the prefixes), and any non-error covariates.
#' @param rep_data  Data frame holding the replicate columns referenced in `formula`.
#' @param family    A binomial GLM family (must be `binomial(link = "logit")`).
#' @param return_details Logical; if TRUE, returns additional internals (xhat, icc, etc.).
#'
#' @return A list with two tidy tables:
#' \itemize{
#'   \item \code{uncorrected}: naive logistic regression estimates.
#'   \item \code{corrected}: sandwich-corrected RC estimates (your final results).
#' }
#' If \code{return_details = TRUE}, also returns intermediate objects.
#'
#' @examples
#' RC_EX_logistic(
#'   formula = comb ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + weight,
#'   main_data = main, rep_data = ers, family = binomial("logit")
#' )
#' @noRd

RC_EX_logistic = function(formula,
                           main_data,
                           rep_data,
                           link = "logistic",
                           return_details = FALSE) {

  # ---- 0) Validate link ----
  if (!is.character(link) || length(link) != 1) {
    stop("`link` must be a single character string such as 'logistic', 'logit', or 'probit'.")
  }

  # Accept "logistic" as an alias for "logit"
  if (link == "logistic") link = "logit"

  # Construct the binomial family with the chosen link
  family = binomial(link = link)

  # ---- 1) Normalize formula to string and parse ----
  formula_str = if (inherits(formula, "formula")) {
    paste(deparse(formula), collapse = "")
  } else {
    as.character(formula)
  }

  parsed = parse_and_extract(
    main_data    = main_data,
    external_data = rep_data,
    formula_str  = formula_str
  )
  # parsed returns: main_data (unchanged), Y, W, z.main, z.rep, r

  # ---- 2) Prepare & standardize ----
  prep = prepare_data_ex(
    z.main = parsed$z.main,
    r      = parsed$r,
    z.rep  = parsed$z.rep,
    W      = parsed$W,
    Y      = parsed$Y
  )
  # prep returns: z.main.std, z.rep.std, r, (optional) W.main.std, Y, indicator,
  #               means = list(z=..., w=...), sds = list(z=..., w=...)

  # ---- 3) Naive logistic regression (for reference) ----
  naive = naive_analysis_ex_log(
    z.main.std  = prep$z.main.std,
    W.main.std  = prep$W.main.std,
    Y           = prep$Y,
    sdz         = prep$sds[["z"]],
    sdw         = prep$sds[["w"]]
  )

  # ---- 4) Regression calibration (xhat, etc.) ----
  rc = reg_calibration_ex_log(
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


  sand <- sandwich_estimator_ex_log(
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
    muz         = muz_std,                     # pass standardized main-study means
    muw         = muw_std,
    fit2        = rc$fit2,
    zbar        = rc$zbar
  )

  out <- list(
    uncorrected = naive[["Naive estimates"]],
    corrected   = sand[["Sandwich Corrected estimates"]]
  )

  if (return_details) {
    out$details <- list(
      parsed      = parsed,
      prepared    = prep,
      rc          = rc
    )
  }

  class(out) <- c("RC_EX_logistic_result", class(out))
  out
}





