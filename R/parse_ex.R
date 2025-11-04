#' Parse Measurement-Error Formula and Extract Data
#'
#' \code{parse_and_extract()} parses a measurement-error formula of the form
#' \code{"Y ~ x(x1, x2) + z(z1, z2, z3) + age + sex"} and extracts the
#' outcome, covariates, error-prone main measurements, replicate measurements,
#' and replicate counts for downstream regression calibration.
#'
#' @param main_data A \code{data.frame} containing the main study data.
#'   Must include the outcome and all non-error-prone variables, and
#'   a single representative column for each measurement-error variable.
#' @param external_data A \code{data.frame} containing the external
#'   reliability/replicate measurements. Must contain all replicate
#'   columns referenced inside the parentheses of \code{formula_str}.
#' @param formula_str A character string specifying the model formula,
#'   e.g. \code{"Y ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + age + sex"}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{main_data}}{The original \code{main_data}.}
#'   \item{\code{Y}}{Outcome vector extracted from \code{main_data}.}
#'   \item{\code{W}}{Matrix of error-free covariates (may be \code{NULL}).}
#'   \item{\code{z.main}}{Matrix of error-prone main measurements, one column per prefix (may be \code{NULL}).}
#'   \item{\code{z.rep}}{Named list of matrices of replicate measurements, one entry per prefix (may be empty).}
#'   \item{\code{r}}{Integer vector giving the number of non-missing replicates for each subject (0 if none).}
#' }
#'
#' @examples
#' # Example main and external data
#' set.seed(1)
#' main_df = data.frame(
#'   Y    = rbinom(5, 1, 0.5),
#'   sbp  = rnorm(5, 120, 15),
#'   age  = rnorm(5, 50, 10),
#'   sex  = sample(c(0, 1), 5, TRUE)
#' )
#' ext_df <- data.frame(
#'   sbp2 = rnorm(5, 120, 15),
#'   sbp3 = rnorm(5, 120, 15)
#' )
#'
#' # Parse formula and extract data
#' parsed <- parse_and_extract(
#'   main_data    = main_df,
#'   external_data = ext_df,
#'   formula_str   = "Y ~ sbp(sbp2, sbp3) + age + sex"
#' )
#'
#' str(parsed)
#'
#' @noRd



parse_and_extract = function(
    main_data,
    external_data,
    formula_str
) {


  parts = strsplit(formula_str, "~")[[1]]
  if (length(parts) != 2) {
    stop("Formula must contain exactly one '~'.")
  }

  lhs = trimws(parts[1])  # outcome name
  rhs = trimws(parts[2])  # everything to the right of '~'


  rhs_terms = strsplit(rhs, "\\+")[[1]] # Splits RHS into individual terms using +
  rhs_terms = trimws(rhs_terms)  # remove extra whitespace


  me_pattern = "^(\\w+)\\s*\\(([^)]*)\\)$" # Identify Measurement-Error Terms
  me_blocks    = list()
  noerr_vars   = character(0)


  for (term in rhs_terms) {
    match = regexpr(me_pattern, term)

    if (match != -1) {
      # e.g., term = "bmi1(bmi1_1, bmi1_2)"
      found  = regmatches(term, match)
      prefix = sub(me_pattern, "\\1", found)    # e.g. "bmi1"
      inside = sub(me_pattern, "\\2", found)    # e.g. "bmi1_1, bmi1_2"

      replicate_cols = strsplit(inside, "\\s*,\\s*")[[1]]
      me_blocks[[prefix]] = replicate_cols

    } else {
      noerr_vars = c(noerr_vars, term) # non-error prone variabel
    }
  }


  if (!lhs %in% names(main_data)) {
    stop(sprintf("Outcome '%s' not found in main_data.", lhs))
  }

  Y = main_data[[lhs]]


  W = NULL

  if (length(noerr_vars) > 0) {
    missing_noerr = setdiff(noerr_vars, names(main_data)) # Checks if any of the variables listed in noerr_vars are missing from main_data
    if (length(missing_noerr) > 0) {
      stop(sprintf(
        "Missing no-error variables in main_data: %s",
        paste(missing_noerr, collapse = ", ")
      ))
    }
    W = as.matrix(main_data[, noerr_vars, drop = FALSE])
  }


  z.main = NULL
  z.rep  = list()

  if (length(me_blocks) > 0) {
    # Build z.main matrix from all ME prefixes
    z_main_list = lapply(
      names(me_blocks),
      function(prefix) {
        if (!prefix %in% names(main_data)) {
          stop(sprintf(
            "Measurement-error variable '%s' not found in main_data.", prefix
          ))
        }
        main_data[[prefix]]
      }
    )

    z.main = do.call(cbind, z_main_list)
    colnames(z.main) = names(me_blocks)

    # Build z.rep list for each ME prefix
    for (prefix in names(me_blocks)) {
      replicate_cols = me_blocks[[ prefix ]]
      missing_reps   = setdiff(replicate_cols, names(external_data))
      if (length(missing_reps) > 0) {
        stop(sprintf(
          "For prefix '%s', replicate columns not found in external_data: %s",
          prefix, paste(missing_reps, collapse=", ")
        ))
      }
      z.rep[[ prefix ]] = as.matrix(external_data[, replicate_cols, drop = FALSE])
    }
  }


  if (length(z.rep) > 0) {
    first_prefix = names(z.rep)[1]
    r = apply(z.rep[[ first_prefix ]], 1, function(x) sum(!is.na(x)))
  } else {
    r = rep(0, nrow(main_data))  # or NA if no ME variables exist
  }


  return(list(
    main_data = main_data,
    Y         = Y,
    W         = W,
    z.main    = z.main,
    z.rep     = z.rep,
    r         = r
  ))
}
