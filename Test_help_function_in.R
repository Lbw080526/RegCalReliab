############################ Parse_in ############################

set.seed(1)
main_df <- data.frame(
  Y     = rbinom(5, 1, 0.5),
  sbp1  = rnorm(5, 120, 15),
  sbp2  = rnorm(5, 120, 15),
  age   = rnorm(5, 50, 10),
  sex   = sample(c(0, 1), 5, TRUE)
)

# Parse formula and extract data (1Z1W example)
parsed <- parse_and_extract_internal(
  main_data   = main_df,
  formula_str = "Y ~ mysbp(sbp1, sbp2) + age + sex"
)

############################ prepare_data_in ############################
set.seed(123)
# Internal study: 5 subjects, 1 exposure with 2 replicates
z <- list(
  sbp = cbind(rnorm(5, 120, 15), rnorm(5, 120, 15))
)

# Replicate counts (each subject has 2 replicates)
r <- rep(2, 5)

# Outcome and optional covariate
Y <- rbinom(5, 1, 0.5)
W <- matrix(rnorm(5), ncol = 1)
colnames(W) <- "age"

# Prepare and standardize data
prep <- prepare_data_in(r = r, z = z, W = W, Y = Y)

############################ naive_analysis_in_log ############################

set.seed(1)
# Simulated internal data: 100 subjects, 2 replicates of 1 exposure
z.rep <- cbind(rnorm(100), rnorm(100))
zbar <- rowMeans(z.rep)
Y <- rbinom(100, 1, plogis(0.5 * zbar))

# Standardize and get SDs
zbar.std <- scale(zbar)
sdz <- sd(zbar)

# Run naive logistic regression ignoring measurement error
res <- naive_analysis_in_log(
  zbar = zbar.std,
  W.std = NULL,
  Y = Y,
  sdz = sdz,
  sdw = NULL
)

############################ reg_calibration_in_log ############################
set.seed(123)
# Internal reliability study: 60 subjects, 2 replicates of 1 exposure
z.rep <- cbind(rnorm(60), rnorm(60))
zbar <- rowMeans(z.rep)
Y <- rbinom(60, 1, plogis(0.4 * zbar))

# Standardize data
zbar.std <- scale(zbar)
sdz <- sd(zbar)
z.std <- list(sbp = scale(z.rep))
r <- rep(2, 60) # each subject has 2 replicates

# Naive covariance (for dimension labels)
naive <- naive_analysis_in_log(Y = Y, zbar = zbar.std, W.std = NULL,
                               sdz = sdz, sdw = NULL)

# Apply regression calibration
fit <- reg_calibration_in_log(
  Y = Y,
  zbar = as.matrix(zbar.std),
  z.std = z.std,
  W.std = NULL,
  muz = mean(zbar),
  muw = NULL,
  sdz = sdz,
  sdw = NULL,
  r = r,
  var1 = naive$var1
)
str(fit)


############################ naive_analysis_in_linear ############################
set.seed(1)
# Simulated internal data: 80 subjects, 2 replicates of 1 exposure
z.rep <- cbind(rnorm(80), rnorm(80))
zbar <- rowMeans(z.rep)
Y <- 2 + 0.5 * zbar + rnorm(80)

# Standardize and get SDs
zbar.std <- scale(zbar)
sdz <- sd(zbar)

# Run naive linear regression ignoring measurement error
res <- naive_analysis_in_linear(
  zbar = zbar.std,
  W.std = NULL,
  Y = Y,
  sdz = sdz,
  sdw = NULL
)


############################ reg_calibration_in_linear ####################
set.seed(123)
# Internal reliability study: 60 subjects, 2 replicates of 1 exposure
z.rep <- cbind(rnorm(60), rnorm(60))
zbar <- rowMeans(z.rep)
Y <- 2 + 0.5 * zbar + rnorm(60)

# Standardize data
zbar.std <- scale(zbar)
sdz <- sd(zbar)
z.std <- list(sbp = scale(z.rep))
r <- rep(2, 60) # each subject has 2 replicates

# Naive covariance (for dimension labels)
naive <- naive_analysis_in_linear(Y = Y, zbar = zbar.std, W.std = NULL,
                                  sdz = sdz, sdw = NULL)

# Apply regression calibration
fit <- reg_calibration_in_linear(
  Y = Y,
  zbar = as.matrix(zbar.std),
  z.std = z.std,
  W.std = NULL,
  muz = mean(zbar),
  muw = NULL,
  sdz = sdz,
  sdw = NULL,
  r = r,
  var1 = naive$var1
)

############################ naive_analysis_in_poisson ####################

set.seed(123)
# Simulated replicate data: 100 subjects, 1 exposure with 2 replicates
z.rep <- cbind(rnorm(100), rnorm(100))
zbar <- rowMeans(z.rep)
Y <- rpois(100, exp(0.3 + 0.5 * zbar))
sdz <- sd(zbar)

# Run naive Poisson regression
res <- naive_analysis_in_poisson(
  Y = Y,
  zbar = scale(zbar),
  W.std = NULL,
  sdz = sdz,
  sdw = NULL
)

############################ reg_calibration_in_poisson ####################

set.seed(123)
# Internal reliability study: 60 subjects, 2 replicates of 1 exposure
z.rep <- cbind(rnorm(60), rnorm(60))
zbar <- rowMeans(z.rep)
Y <- rpois(60, lambda = exp(0.2 + 0.5 * zbar))

# Standardize data
zbar.std <- scale(zbar)
sdz <- sd(zbar)
z.std <- list(sbp = scale(z.rep))
r <- rep(2, 60)

# Naive covariance (for dimension labels)
naive <- naive_analysis_in_poisson(Y = Y, zbar = zbar.std,
                                   W.std = NULL, sdz = sdz, sdw = NULL)

# Apply regression calibration
fit <- reg_calibration_in_poisson(
  Y = Y,
  zbar = as.matrix(zbar.std),
  z.std = z.std,
  W.std = NULL,
  muz = mean(zbar),
  muw = NULL,
  sdz = sdz,
  sdw = NULL,
  r = r,
  var1 = naive$var1
)





set.seed(123)
add_err <- function(v, sd = sqrt(0.4)) v + rnorm(length(v), 0, sd)

## --- Example 1: Internal 1Z 0W ---
x <- rnorm(3000)
z <- rbind(
  cbind(add_err(x[1:1500]), NA, NA, NA),
  cbind(add_err(x[1501:2000]), add_err(x[1501:2000]), NA, NA),
  cbind(add_err(x[2001:2400]), add_err(x[2001:2400]), add_err(x[2001:2400]), NA),
  cbind(add_err(x[2401:3000]), add_err(x[2401:3000]), add_err(x[2401:3000]), add_err(x[2401:3000]))
)
colnames(z) <- paste0("z_", 1:4)
Y <- rbinom(3000, 1, plogis(-2.65 + log(1.5) * x))
main_data <- data.frame(Y, z)
res1 <- RC_InReliab(Y ~ myz(z_1, z_2, z_3, z_4),
                    main_data = main_data,
                    link = "logistic")
res1$corrected

## --- Example 2: Internal 1Z 1W ---
x  <- rnorm(3000)
W1 <- rnorm(3000)
z <- rbind(
  cbind(add_err(x[1:1500]), NA, NA, NA),
  cbind(add_err(x[1501:2000]), add_err(x[1501:2000]), NA, NA),
  cbind(add_err(x[2001:2400]), add_err(x[2001:2400]), add_err(x[2001:2400]), NA),
  cbind(add_err(x[2401:3000]), add_err(x[2401:3000]), add_err(x[2401:3000]), add_err(x[2401:3000]))
)
colnames(z) <- paste0("z_", 1:4)
Y <- rbinom(3000, 1, plogis(-2.65 + log(1.5) * x + 0.5 * W1))
main_data <- data.frame(Y, z, W1)
res2 <- RC_InReliab(Y ~ myz(z_1, z_2, z_3, z_4) + W1,
                    main_data = main_data,
                    link = "logistic")
res2$corrected

## --- Example 3: Internal 2Z 0W ---
x <- mgcv::rmvn(3000, c(0,0), matrix(c(1,0.3,0.3,1), 2))
z1 <- rbind(
  cbind(add_err(x[1:1500, 1]), NA, NA, NA),
  cbind(add_err(x[1501:2000, 1]), add_err(x[1501:2000, 1]), NA, NA),
  cbind(add_err(x[2001:2400, 1]), add_err(x[2001:2400, 1]), add_err(x[2001:2400, 1]), NA),
  cbind(add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]), add_err(x[2401:3000, 1]))
)
colnames(z1) <- paste0("z1_", 1:4)
z2 <- rbind(
  cbind(add_err(x[1:1500, 2]), NA, NA, NA),
  cbind(add_err(x[1501:2000, 2]), add_err(x[1501:2000, 2]), NA, NA),
  cbind(add_err(x[2001:2400, 2]), add_err(x[2001:2400, 2]), add_err(x[2001:2400, 2]), NA),
  cbind(add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]), add_err(x[2401:3000, 2]))
)
colnames(z2) <- paste0("z2_", 1:4)
Y <- rbinom(3000, 1, plogis(-2.65 + log(1.5) * rowSums(x)))
main_data <- data.frame(Y, z1, z2)
res3 <- RC_InReliab(
  Y ~ myz1(z1_1, z1_2, z1_3, z1_4) + myz2(z2_1, z2_2, z2_3, z2_4),
  main_data = main_data,
  link = "logistic")
res3$corrected
