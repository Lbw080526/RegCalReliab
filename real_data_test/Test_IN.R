library(Regcal)
library(mgcv)

# Add error function
add_err <- function(v, sd = sqrt(0.4)) v + rnorm(length(v), 0, sd)

# One replicate
simulate_once_internal <- function() {
  x <- mgcv::rmvn(3000, c(0, 0), matrix(c(1, 0.3, 0.3, 1), 2))
  r <- c(rep(1, 1500), rep(2, 500), rep(3, 400), rep(4, 600))

  z <- list(
    rbind(
      cbind(add_err(x[1:1500, 1]), NA, NA, NA),
      cbind(add_err(x[1501:2000, 1]), add_err(x[1501:2000, 1]), NA, NA),
      cbind(add_err(x[2001:2400, 1]), add_err(x[2001:2400, 1]), add_err(x[2001:2400, 1]), NA),
      cbind(add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]))
    ),
    rbind(
      cbind(add_err(x[1:1500, 2]), NA, NA, NA),
      cbind(add_err(x[1501:2000, 2]), add_err(x[1501:2000, 2]), NA, NA),
      cbind(add_err(x[2001:2400, 2]), add_err(x[2001:2400, 2]), add_err(x[2001:2400, 2]), NA),
      cbind(add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]))
    )
  )

  p <- exp(-2.3 + log(1.5) * rowSums(x)) / (1 + exp(-2.3 + log(1.5) * rowSums(x)))
  Y <- rbinom(3000, 1, p)

  results <- RC_IN(model = "logistic", r = r, z = z, W = NULL, Y = Y)
  return(results)
}

# 让yu看看名字, suggest: RegcalRELI
# Run simulation
Nrep <- 10
results_list <- replicate(Nrep, simulate_once_internal(), simplify = FALSE)

# Rename rows
rename_naive_rows <- function(tab) {
  non_int <- setdiff(rownames(tab), "(Intercept)")
  if (length(non_int) == 2)
    rownames(tab)[match(non_int, rownames(tab))] <- c("z1_df", "z2_df")
  tab
}
rename_hat_rows <- function(tab) {
  hats <- grep("_hat$", rownames(tab))
  if (length(hats) == 2)
    rownames(tab)[hats] <- c("z1_hat", "z2_hat")
  tab
}

results_list <- lapply(results_list, function(res) {
  res$uncorrected <- rename_naive_rows(res$uncorrected)
  res$corrected   <- rename_hat_rows(res$corrected)
  res
})

# Extract rows
get_row <- function(tab, rowname) tab[rowname, , drop = FALSE]

# Combine rows
naive_z1     <- do.call(rbind, lapply(results_list, \(res) get_row(res$uncorrected, "z1_df")))
corrected_z1 <- do.call(rbind, lapply(results_list, \(res) get_row(res$corrected,   "z1_hat")))
naive_z2     <- do.call(rbind, lapply(results_list, \(res) get_row(res$uncorrected, "z2_df")))
corrected_z2 <- do.call(rbind, lapply(results_list, \(res) get_row(res$corrected,   "z2_hat")))

# Averages
avg_naive_z1     <- colMeans(naive_z1)
avg_corrected_z1 <- colMeans(corrected_z1)
avg_naive_z2     <- colMeans(naive_z2)
avg_corrected_z2 <- colMeans(corrected_z2)

cat("\nAverage across", Nrep, "replicates (intercept omitted):\n\n")
print(rbind(
  Naive_Z1     = avg_naive_z1,
  Corrected_Z1 = avg_corrected_z1
))
print(rbind(
  Naive_Z2     = avg_naive_z2,
  Corrected_Z2 = avg_corrected_z2
))

# --- Coverage
OR_true <- 1.5

inside_ci <- function(tab, rowname, truth = OR_true) {
  if (rowname %in% rownames(tab)) {
    ci <- tab[rowname, c("CI.low", "CI.high")]
    return(ci[1] <= truth && truth <= ci[2])
  }
  return(NA)
}

coverage <- function(component, rowname) {
  mean(sapply(results_list, \(res) inside_ci(res[[component]], rowname)), na.rm = TRUE) * 100
}

# Z1
cov_z1_naive = coverage("uncorrected", "z1_df")
cov_z1_corr  = coverage("corrected",   "z1_hat")

# Z2
cov_z2_naive = coverage("uncorrected", "z2_df")
cov_z2_corr  = coverage("corrected",   "z2_hat")

cat("\nCoverage of TRUE OR = 1.5 for error-prone exposure z1:\n")
cat(sprintf("  • Naive     : %5.1f %%\n", cov_z1_naive))
cat(sprintf("  • Corrected : %5.1f %%\n\n", cov_z1_corr))

cat("Coverage of TRUE OR = 1.5 for error-prone exposure z2:\n")
cat(sprintf("  • Naive     : %5.1f %%\n", cov_z2_naive))
cat(sprintf("  • Corrected : %5.1f %%\n", cov_z2_corr))
