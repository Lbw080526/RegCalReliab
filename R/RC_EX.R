#' Unified Regression Calibration Wrapper (External Reliability Study)
#'
#' @description
#' A single formula interface for regression calibration in external reliability
#' studies. The user simply specifies `link = "linear"`, `"logistic"`, or `"log"`,
#' and the wrapper selects the appropriate model:
#'   * `"linear"`   → Gaussian (identity link)
#'   * `"logistic"` → Binomial (logit link)
#'   * `"log"`      → Poisson (log link)
#'
#' @param formula A formula or character string such as
#'   "Y ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + weight".
#'   Terms of the form `var(rep1, rep2, ...)` are treated as error-prone exposures
#'   with replicates in `rep_data`; other terms are treated as covariates W.
#' @param main_data Data frame holding the outcome, error-prone exposures, and covariates.
#' @param rep_data  Data frame holding replicate columns referenced in `formula`.
#' @param link Character; one of `"linear"`, `"logistic"`, or `"log"`.
#' @param return_details Logical; if TRUE, return parsed, prepared, and RC internals.
#'
#' @return A list with
#'   * `uncorrected`: naive regression estimates
#'   * `corrected`  : sandwich-corrected regression calibration estimates
#'   plus optional `details` if `return_details = TRUE`.
#'
#' @export
RC_ExReliab <- function(formula,
                        main_data,
                        rep_data,
                        link = c("linear", "logistic", "log"),
                        return_details = FALSE) {

  # ---- 0) Validate link ----
  link <- match.arg(link)
  # Normalize link so each sub-wrapper can rely on its expected input
  if (link == "linear")   family_name <- "gaussian"
  if (link == "logistic") family_name <- "binomial"
  if (link == "log")      family_name <- "poisson"

  # ---- 1) Dispatch to the appropriate model wrapper ----
  if (link == "linear") {
    out <- RC_EX_Linear(
      formula        = formula,
      main_data      = main_data,
      rep_data       = rep_data,
      link           = "linear",
      return_details = return_details
    )

  } else if (link == "logistic") {
    out <- RC_EX_logistic(
      formula        = formula,
      main_data      = main_data,
      rep_data       = rep_data,
      link           = "logistic",
      return_details = return_details
    )

  } else if (link == "log") {
    out <- RC_EX_Poisson(
      formula        = formula,
      main_data      = main_data,
      rep_data       = rep_data,
      link           = "log",
      return_details = return_details
    )
  }

  out
}
