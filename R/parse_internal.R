#' Parse Measurement-Error Formula and Extract Data (Internal Reliability)
#'
#' \code{parse_and_extract_internal()} parses a measurement-error formula of the
#' form \code{"Y ~ x(x1, x2) + z(z1, z2, z3) + age + sex"} for \emph{internal}
#' reliability study designs. It extracts the outcome, covariates, error-prone
#' replicate measurements, and replicate counts for downstream regression
#' calibration.
#'
#' @param main_data A \code{data.frame} containing the main study data.
#'   Must include the outcome, all non-error-prone variables, and all replicate
#'   columns referenced inside the parentheses of \code{formula_str}.
#' @param formula_str A character string specifying the model formula,
#'   e.g. \code{"Y ~ sbp(sbp1, sbp2) + chol(chol1, chol2, chol3) + age + sex"}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{Y}}{Outcome vector extracted from \code{main_data}.}
#'   \item{\code{W}}{Matrix of error-free covariates (may be \code{NULL}).}
#'   \item{\code{z}}{Named list of matrices of replicate measurements for each
#'     error-prone exposure (may be empty).}
#'   \item{\code{r}}{Integer vector giving the number of non-missing replicates
#'     for each subject (0 if none).}
#' }
#'
#' @examples
#' # Example main data with replicates
#' set.seed(1)
#' main_df <- data.frame(
#'   Y     = rbinom(5, 1, 0.5),
#'   sbp1  = rnorm(5, 120, 15),
#'   sbp2  = rnorm(5, 120, 15),
#'   age   = rnorm(5, 50, 10),
#'   sex   = sample(c(0, 1), 5, TRUE)
#' )
#'
#' # Parse formula and extract data (1Z1W example)
#' parsed <- parse_and_extract_internal(
#'   main_data   = main_df,
#'   formula_str = "Y ~ mysbp(sbp1, sbp2) + age + sex"
#' )
#'
#' str(parsed)
#'
#' @noRd



parse_and_extract_internal = function(main_data, formula_str) {

  parts = strsplit(formula_str, "~")[[1]]
  if (length(parts) != 2) {
    stop("Formula must contain exactly one '~'.")
  }
  lhs = trimws(parts[1])
  rhs = trimws(parts[2])

  if (!lhs %in% names(main_data)) {
    stop(sprintf("Outcome '%s' not found in main_data.", lhs))
  }
  Y = main_data[[lhs]]


  rhs_terms = strsplit(rhs, "\\+")[[1]]
  rhs_terms = trimws(rhs_terms)

  noerr_vars = character(0)
  me_blocks  = list()


  re_prefix_paren = "^([A-Za-z]\\w*)\\((.*)\\)$"

  for (term in rhs_terms) {

    term = trimws(term)

    if (grepl(re_prefix_paren, term)) {
      prefix = sub(re_prefix_paren, "\\1", term)  # e.g. "BMI"
      inside = sub(re_prefix_paren, "\\2", term)  # e.g. "BMI_1, BMI_2, BMI_3"

      replicate_cols = strsplit(inside, "\\s*,\\s*")[[1]]

      me_blocks[[prefix]] = replicate_cols

    } else if (grepl("^\\(.*\\)$", term)) {
      inside = sub("^\\((.*)\\)$", "\\1", term)
      replicate_cols = strsplit(inside, "\\s*,\\s*")[[1]]

      prefix = sub("^([A-Za-z]+)_.*", "\\1", replicate_cols[1])

      me_blocks[[prefix]] = replicate_cols

    } else {
      noerr_vars = c(noerr_vars, term)
    }
  }

  W = NULL
  if (length(noerr_vars) > 0) {
    missing_noerr = setdiff(noerr_vars, names(main_data))
    if (length(missing_noerr) > 0) {
      stop(sprintf(
        "Missing no-error variables in main_data: %s",
        paste(missing_noerr, collapse = ", ")
      ))
    }
    W = as.matrix(main_data[, noerr_vars, drop = FALSE])
  }

  z = list()
  if (length(me_blocks) > 0) {
    for (prefix in names(me_blocks)) {

      replicate_cols = me_blocks[[prefix]]
      # Check that replicate columns actually exist
      missing_reps = setdiff(replicate_cols, names(main_data))
      if (length(missing_reps) > 0) {
        stop(sprintf(
          "For prefix '%s', replicate columns not found in main_data: %s",
          prefix, paste(missing_reps, collapse = ", ")
        ))
      }

      z[[prefix]] = as.matrix(main_data[, replicate_cols, drop = FALSE])
    }

    first_prefix = names(z)[1]
    r = apply(z[[first_prefix]], 1, function(x) sum(!is.na(x)))

  } else {
    r = rep(0, nrow(main_data))
  }

  return(list(
    Y = Y,
    W = W,
    z = z,
    r = r
  ))
}
