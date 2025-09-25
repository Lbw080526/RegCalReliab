############################ Parse_ex ############################
set.seed(1)
main_df <- data.frame(
  Y    = rbinom(5, 1, 0.5),
  sbp  = rnorm(5, 120, 15),
  age  = rnorm(5, 50, 10),
  sex  = sample(c(0, 1), 5, TRUE)
)
ext_df <- data.frame(
  sbp2 = rnorm(5, 120, 15),
  sbp3 = rnorm(5, 120, 15)
)

# Parse formula and extract data
parsed <- parse_and_extract(
  main_data    = main_df,
  external_data = ext_df,
  formula_str   = "Y ~ sbp(sbp2, sbp3) + age + sex"
)


############################ Prepare_EX ############################
set.seed(123)
# Main study: 6 subjects, 2 error-prone exposures
z.main <- matrix(rnorm(6 * 2), nrow = 6, ncol = 2)
colnames(z.main) <- c("sbp", "chol")

# Reliability study: 4 subjects, each with 2 replicates
z.rep <- list(
  sbp  = matrix(rnorm(4 * 2), nrow = 4),
  chol = matrix(rnorm(4 * 2), nrow = 4)
)

# Replicate counts
r <- rep(2, 4)

# Optional covariates and binary outcome
W <- matrix(rnorm(6 * 2), nrow = 6)
colnames(W) <- c("age", "sex")
Y <- rbinom(6, 1, 0.5)

# Prepare and standardize data
prep <- prepare_data_ex(z.main = z.main,
                        r = r,
                        z.rep = z.rep,
                        W = W,
                        Y = Y)


############################ naive_ex_log ############################
set.seed(1)
# Simulated main-study data: 100 subjects, 1 exposure
z <- matrix(rnorm(100), ncol = 1)
colnames(z) <- "sbp"
Y <- rbinom(100, 1, plogis(0.3 * z))
sdz <- apply(z, 2, sd)

# Run naive logistic regression ignoring measurement error
res <- naive_analysis_ex_log(
  z.main.std = scale(z),
  W.main.std = NULL,
  Y = Y,
  sdz = sdz,
  sdw = NULL
)

############################ reg_calibration_ex_log ############################

set.seed(123)
# Simulated main-study data: 80 subjects, 1 exposure
z.main <- matrix(rnorm(80), ncol = 1)
colnames(z.main) <- "sbp"
Y <- rbinom(80, 1, plogis(0.3 * z.main))

# Reliability study: 40 subjects, 2 replicates
z.rep <- list(sbp = matrix(rnorm(40 * 2), nrow = 40))
r <- c(rep(1, 80), rep(2, 40)) # replicate counts
indicator <- c(rep(1, 80), rep(0, 40))

# Standardize data
sdz <- apply(z.main, 2, sd)
z.main.std <- scale(z.main)
z.rep.std <- list(sbp = scale(z.rep$sbp))

# Apply regression calibration
fit <- reg_calibration_ex_log(
  z.main.std = z.main.std,
  z.rep.std  = z.rep.std,
  r          = r,
  W.main.std = NULL,
  Y          = Y,
  muz        = colMeans(z.main),
  muw        = NULL,
  sdz        = sdz,
  sdw        = NULL,
  indicator  = indicator
)

########################## naive_analysis_ex_linear #########################
set.seed(1)
# Simulated main-study data: 100 subjects, 1 exposure
z <- matrix(rnorm(100), ncol = 1)
colnames(z) <- "sbp"
Y <- 2 + 0.5 * z + rnorm(100)
sdz <- apply(z, 2, sd)

# Run naive linear regression ignoring measurement error
res <- naive_analysis_ex_linear(
  z.main.std = scale(z),
  W.main.std = NULL,
  Y = Y,
  sdz = sdz,
  sdw = NULL
)


########################## reg_calibration_ex_linear #########################
set.seed(1)
# Simulated main-study data: 80 subjects, 1 exposure
z <- matrix(rnorm(80), ncol = 1)
colnames(z) <- "sbp"
Y <- 2 + 0.5 * z + rnorm(80)

# Reliability study: 40 subjects, 2 replicates
z.rep <- list(sbp = matrix(rnorm(40 * 2), nrow = 40))
r <- c(rep(1, 80), rep(2, 40))
indicator <- c(rep(1, 80), rep(0, 40))

# Standardize data
sdz <- apply(z, 2, sd)
z.main.std <- scale(z)
z.rep.std <- list(sbp = scale(z.rep$sbp))

# Run regression calibration
fit <- reg_calibration_ex_linear(
  z.main.std = z.main.std,
  z.rep.std  = z.rep.std,
  r          = r,
  W.main.std = NULL,
  Y          = Y,
  muz        = colMeans(z),
  muw        = NULL,
  sdz        = sdz,
  sdw        = NULL,
  indicator  = indicator
)



########################## naive_analysis_ex_poisson #########################

set.seed(1)
# Simulated main-study data: 100 subjects, 1 exposure
z <- matrix(rnorm(100), ncol = 1)
colnames(z) <- "exposure"
Y <- rpois(100, lambda = exp(0.2 + 0.3 * z))
sdz <- apply(z, 2, sd)

# Run naive Poisson regression ignoring measurement error
res <- naive_analysis_ex_poisson(
  z.main.std = scale(z),
  W.main.std = NULL,
  Y = Y,
  sdz = sdz,
  sdw = NULL
)


########################## naive_analysis_ex_poisson #########################
set.seed(123)
# Simulated main-study data: 80 subjects, 1 exposure
z.main <- matrix(rnorm(80), ncol = 1)
colnames(z.main) <- "sbp"
Y <- rpois(80, lambda = exp(0.3 * z.main))

# Reliability study: 40 subjects, 2 replicates
z.rep <- list(sbp = matrix(rnorm(40 * 2), nrow = 40))
r <- c(rep(1, 80), rep(2, 40)) # replicate counts
indicator <- c(rep(1, 80), rep(0, 40))

# Standardize data
sdz <- apply(z.main, 2, sd)
z.main.std <- scale(z.main)
z.rep.std <- list(sbp = scale(z.rep$sbp))

# Apply regression calibration for Poisson regression
fit <- reg_calibration_ex_poisson(
  z.main.std = z.main.std,
  z.rep.std  = z.rep.std,
  r          = r,
  W.main.std = NULL,
  Y          = Y,
  muz        = colMeans(z.main),
  muw        = NULL,
  sdz        = sdz,
  sdw        = NULL,
  indicator  = indicator
)


########################## RC_ExReliab #########################
x <- rnorm(3000)
z.main <- x[1:1500] + rnorm(1500, 0, sqrt(0.4))
z_rep <- rbind(
  cbind(add_err(x[1501:2000]), add_err(x[1501:2000]), NA, NA),
  cbind(add_err(x[2001:2400]), add_err(x[2001:2400]), add_err(x[2001:2400]), NA),
  cbind(add_err(x[2401:3000]), add_err(x[2401:3000]), add_err(x[2401:3000]), add_err(x[2401:3000]))
)
colnames(z_rep) <- paste0("z_", 1:4)
Y <- rbinom(1500, 1, plogis(-2.3 + log(1.5) * x[1:1500]))
main_data <- data.frame(Y = Y, z = z.main)
rep_data  <- data.frame(z_rep, check.names = FALSE)
res1 <- RC_ExReliab(Y ~ z(z_1, z_2, z_3, z_4), main_data, rep_data, link = "logistic")
res1$corrected

x <- rnorm(3000)
W_main <- rnorm(1500)
W_rep  <- rnorm(1500)
z.main <- x[1:1500] + rnorm(1500, 0, sqrt(0.4))
z_rep <- rbind(
  cbind(add_err(x[1501:2000]), add_err(x[1501:2000]), NA, NA),
  cbind(add_err(x[2001:2400]), add_err(x[2001:2400]), add_err(x[2001:2400]), NA),
  cbind(add_err(x[2401:3000]), add_err(x[2401:3000]), add_err(x[2401:3000]), add_err(x[2401:3000]))
)
colnames(z_rep) <- paste0("z_", 1:4)
Y <- rbinom(1500, 1, plogis(-2.3 + log(1.5) * x[1:1500] + 0.5 * W_main))
main_data <- data.frame(Y = Y, z = z.main, W = W_main)
rep_data  <- data.frame(z_rep, W = W_rep, check.names = FALSE)
res2 <- RC_ExReliab(Y ~ z(z_1, z_2, z_3, z_4) + W, main_data, rep_data, link = "logistic")
res2$corrected

x <- mgcv::rmvn(3000, c(0, 0), matrix(c(1, 0.3, 0.3, 1), 2))
z.main <- x[1:1500, ] + matrix(rnorm(1500 * 2, 0, sqrt(0.4)), 1500, 2)
colnames(z.main) <- c("z1", "z2")
z1_rep <- rbind(
  cbind(add_err(x[1501:2000, 1]), add_err(x[1501:2000, 1]), NA, NA),
  cbind(add_err(x[2001:2400, 1]), add_err(x[2001:2400, 1]), add_err(x[2001:2400, 1]), NA),
  cbind(add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]))
)
colnames(z1_rep) <- paste0("z1_", 1:4)
z2_rep <- rbind(
  cbind(add_err(x[1501:2000, 2]), add_err(x[1501:2000, 2]), NA, NA),
  cbind(add_err(x[2001:2400, 2]), add_err(x[2001:2400, 2]), add_err(x[2001:2400, 2]), NA),
  cbind(add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]))
)
colnames(z2_rep) <- paste0("z2_", 1:4)
Y <- rbinom(1500, 1, plogis(-2.3 + log(1.5) * rowSums(x[1:1500, ])))
main_data <- data.frame(Y = Y, z1 = z.main[, 1], z2 = z.main[, 2])
rep_data  <- data.frame(z1_rep, z2_rep, check.names = FALSE)
res3 <- RC_ExReliab(
  Y ~ z1(z1_1, z1_2, z1_3, z1_4) + z2(z2_1, z2_2, z2_3, z2_4),
  main_data, rep_data, link = "logistic"
)
res3$corrected





