library(Regcal)

# Function to add error to the data
add_err = function(v, sd = sqrt(0.4)) v + rnorm(length(v), 0, sd)

simulate_once <- function() {
  # 1) Simulate data
  x  = mgcv::rmvn(3000, c(0, 0), matrix(c(1, 0.3, 0.3, 1), 2))
  z.main = x[1:1500, 1:2] + matrix(rnorm(1500 * 2, 0, sqrt(0.4)), 1500, 2)

  r = c(rep(2, 500), rep(3, 400), rep(4, 600))

  z.rep = list(
    rbind(
      cbind(add_err(x[1501:2000, 1]), add_err(x[1501:2000, 1]), NA, NA),
      cbind(add_err(x[2001:2400, 1]), add_err(x[2001:2400, 1]), add_err(x[2001:2400, 1]), NA),
      cbind(add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]))
    ),
    rbind(
      cbind(add_err(x[1501:2000, 2]), add_err(x[1501:2000, 2]), NA, NA),
      cbind(add_err(x[2001:2400, 2]), add_err(x[2001:2400, 2]), add_err(x[2001:2400, 2]), NA),
      cbind(add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]))
    )
  )

  p = exp(-2.3 + log(1.5) * rowSums(x[1:1500,])) / (1 + exp(-2.3 + log(1.5) * rowSums(x[1:1500,])))
  Y = rbinom(1500, 1, p)

  # Run regression calibration for logistic model
  results = RC_EX(
    model = "logistic",
    z.main = z.main,
    r = r,
    z.rep = z.rep,
    W = NULL,
    Y = Y,
    muz = colMeans(z.main),
    muw = NULL,
    sdz = apply(z.main, 2, sd),
    sdw = NULL,
    indicator = c(rep(1, length(Y)), rep(0, length(r)))
  )

  return(results)
}

# --- Run simulation
set.seed(2025)
Nrep = 1000
results_list = replicate(Nrep, simulate_once(), simplify = FALSE)

# --- Row name helpers
rename_hat_rows = function(tab) {
  hats = grep("_hat$", rownames(tab))
  if (length(hats) == 2)
    rownames(tab)[hats] = c("z1_hat", "z2_hat")
  tab
}
rename_naive_rows = function(tab) {
  non_int = setdiff(rownames(tab), "(Intercept)")
  if (length(non_int) == 2)
    rownames(tab)[match(non_int, rownames(tab))] = c("z1_df", "z2_df")
  tab
}

# --- Apply row renaming
results_list = lapply(results_list, function(res) {
  res$uncorrected = rename_naive_rows(res$uncorrected)
  res$corrected   = rename_hat_rows(res$corrected)
  res
})

# --- Helper
get_row = function(tab, rowname) tab[rowname, , drop = FALSE]

# --- Aggregate results
naive_z1     = do.call(rbind, lapply(results_list, function(res) get_row(res$uncorrected, "z1_df")))
corrected_z1 = do.call(rbind, lapply(results_list, function(res) get_row(res$corrected, "z1_hat")))

naive_z2     = do.call(rbind, lapply(results_list, function(res) get_row(res$uncorrected, "z2_df")))
corrected_z2 = do.call(rbind, lapply(results_list, function(res) get_row(res$corrected, "z2_hat")))

avg_naive_z1     = colMeans(naive_z1)
avg_corrected_z1 = colMeans(corrected_z1)

avg_naive_z2     = colMeans(naive_z2)
avg_corrected_z2 = colMeans(corrected_z2)

# --- Print average results
cat("\nAverage across", Nrep, "replicates (intercept omitted):\n\n")
print(rbind(
  Naive_Z1    = avg_naive_z1,
  Corrected_Z1 = avg_corrected_z1
))
print(rbind(
  Naive_Z2    = avg_naive_z2,
  Corrected_Z2 = avg_corrected_z2
))

# --- Coverage
OR_true = 1.5

inside_ci = function(tab, rowname, truth = OR_true) {
  if (rowname %in% rownames(tab)) {
    ci = tab[rowname, c("CI.low", "CI.high")]
    return(ci[1] <= truth && truth <= ci[2])
  }
  return(NA)
}

coverage = function(component, rowname) {
  mean(sapply(results_list, function(res) inside_ci(res[[component]], rowname)), na.rm = TRUE) * 100
}

# Z1
cov_z1_naive = coverage("uncorrected", "z1_df")
cov_z1_corr  = coverage("corrected",   "z1_hat")

# Z2
cov_z2_naive = coverage("uncorrected", "z2_df")
cov_z2_corr  = coverage("corrected",   "z2_hat")

# --- Print coverage
cat("\nCoverage of TRUE OR = 1.5 for error-prone exposure z1:\n")
cat(sprintf("  • Naive     : %5.1f %%\n", cov_z1_naive))
cat(sprintf("  • Corrected : %5.1f %%\n\n", cov_z1_corr))

cat("Coverage of TRUE OR = 1.5 for non-error covariate z2:\n")
cat(sprintf("  • Naive     : %5.1f %%\n", cov_z2_naive))
cat(sprintf("  • Corrected : %5.1f %%\n", cov_z2_corr))




# For z1 (error-prone exposure)
sd_naive_z1    <- sd(naive_z1[, "Estimate"],      na.rm = TRUE)
sd_corrected_z1 <- sd(corrected_z1[, "Estimate"], na.rm = TRUE)


# For z2 (non-error covariate)
sd_naive_z2    <- sd(naive_z2[, "Estimate"],      na.rm = TRUE)
sd_corrected_z2 <- sd(corrected_z2[, "Estimate"], na.rm = TRUE)


cat("\nStandard Deviations of Estimates Across", Nrep, "Simulations:\n")

cat("\nFor z1 (beta1, error-prone):\n")
cat(sprintf("  Naive:     %8.5f\n", sd_naive_z1))
cat(sprintf("  Corrected: %8.5f\n", sd_corrected_z1))

cat("\nFor z2 (beta2, non-error):\n")
cat(sprintf("  Naive:     %8.5f\n", sd_naive_z2))
cat(sprintf("  Corrected: %8.5f\n", sd_corrected_z2))






