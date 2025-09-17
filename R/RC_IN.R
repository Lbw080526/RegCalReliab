#' Unified Regression Calibration Wrapper (Internal Reliability Study)
#'
#' @description
#' A single formula interface for regression calibration in internal reliability
#' studies. The user simply specifies `link = "linear"`, `"logistic"`, or `"log"`,
#' and the wrapper selects the appropriate model:
#'   * `"linear"`   → Gaussian (identity link)
#'   * `"logistic"` → Binomial (logit link)
#'   * `"log"`      → Poisson (log link)
#'
#' @param formula A formula or character string such as
#'   "Y ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + weight".
#'   Terms of the form `var(rep1, rep2, ...)` are treated as error-prone exposures
#'   with replicates in `main_data`; other terms are treated as covariates W.
#' @param main_data Data frame holding the outcome, replicate error-prone exposures,
#'   and any covariates.
#' @param link Character; one of `"linear"`, `"logistic"`, or `"log"`.
#' @param return_details Logical; if TRUE, return parsed, prepared, and RC internals.
#'
#' @return A list with
#'   * `uncorrected`: naive regression estimates
#'   * `corrected`  : sandwich-corrected regression calibration estimates
#'   plus optional `details` if `return_details = TRUE`.
#'
#' @export

RC_InReliab <- function(formula,
                        main_data,
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
    out <- RC_IN_Linear(
      formula        = formula,
      main_data      = main_data,
      return_details = return_details,
      link = "linear"
    )

  } else if (link == "logistic") {
    out <- RC_IN_Logistic(
      formula        = formula,
      main_data      = main_data,
      return_details = return_details,
      link = "logistic"
    )

  } else if (link == "log") {
    out <- RC_IN_Poisson(
      formula        = formula,
      main_data      = main_data,
      return_details = return_details,
      link = "log"
    )
  }

  out
}
