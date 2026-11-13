library(tidyverse)

# Fast version of fit_agreement() for balanced designs, used in the Monte Carlo studies.
#
# With I patients and J measurements each, the REML estimates of the random-intercept
# model D_ij = beta0 + u_i + e_ij have a closed form (balanced one-way ANOVA):
#   sigma_e^2 = MS_W,  sigma_u^2 = (MS_B - MS_W) / J,  beta0 = grand mean.
# This is much faster than lmer(), which matters for 1,000 datasets x 999 bootstrap refits.
#
# The differences are passed as a matrix with one row per patient and one column per measurement.

# Classical analysis and additive/multiplicative model from the ANOVA sums of squares.
# grand_mean, ss_within and ss_between can be vectors (one value per bootstrap replicate).
agreement_from_ss <- function(grand_mean, ss_within, ss_between, I, J, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  n <- I * J
  ms_within <- ss_within / (I * (J - 1))
  ms_between <- ss_between / (I - 1)

  # Classical analysis: SD of all n differences, as if they were independent
  sd_classical <- sqrt((ss_within + ss_between) / (n - 1))

  # Random-intercept model (REML). Variances cannot be negative: if MS_B < MS_W,
  # the REML optimum is sigma_u^2 = 0, and sigma_e^2 is then re-estimated from all
  # n differences (it is not MS_W any more).
  boundary <- ms_between < ms_within
  var_patient <- if_else(boundary, 0, (ms_between - ms_within) / J)
  var_residual <- if_else(boundary, sd_classical^2, ms_within)
  sd_total <- sqrt(var_patient + var_residual)

  tibble(
    # Random-intercept model, same names as fit_agreement()
    bias = grand_mean,
    se = sqrt(ms_between / n),  # = SD of the I patient means / sqrt(I)
    sd_patient = sqrt(var_patient),
    sd_residual = sqrt(var_residual),
    boundary = boundary,
    lower_loa = grand_mean - z * sd_total,
    upper_loa = grand_mean + z * sd_total,
    # Classical analysis (the bias is the same mean in a balanced design)
    classical_bias = grand_mean,
    classical_se = sd_classical / sqrt(n),
    classical_lower_loa = grand_mean - z * sd_classical,
    classical_upper_loa = grand_mean + z * sd_classical
  )
}

# Fit to one dataset (d: I x J matrix of differences)
fit_balanced <- function(d) {
  I <- nrow(d)
  J <- ncol(d)
  patient_mean <- rowMeans(d)
  grand_mean <- mean(patient_mean)

  agreement_from_ss(
    grand_mean = grand_mean,
    ss_within = sum((d - patient_mean)^2),  # subtracts each patient's own mean
    ss_between = J * sum((patient_mean - grand_mean)^2),
    I = I, J = J
  )
}

# Patient bootstrap: resample whole patients (rows of d) n_boot times.
# The fit only needs each patient's mean and within-patient sum of squares, so these
# are computed once and each bootstrap replicate just resamples I pairs of numbers.
# Returns one row per bootstrap replicate (classical and model estimates).
bootstrap_patients <- function(d, n_boot) {
  I <- nrow(d)
  J <- ncol(d)
  patient_mean <- rowMeans(d)
  patient_ss <- rowSums((d - patient_mean)^2)

  # Column h holds the I patients drawn in bootstrap replicate h
  draws <- matrix(sample.int(I, I * n_boot, replace = TRUE), nrow = I)
  boot_mean <- matrix(patient_mean[draws], nrow = I)
  boot_ss <- matrix(patient_ss[draws], nrow = I)

  grand_mean <- colMeans(boot_mean)
  agreement_from_ss(
    grand_mean = grand_mean,
    ss_within = colSums(boot_ss),
    # sum of (patient mean - grand mean)^2 = sum of squares - I * grand mean^2
    ss_between = J * (colSums(boot_mean^2) - I * grand_mean^2),
    I = I, J = J
  )
}

# Observation bootstrap for the classical analysis: resample single measurements,
# ignoring the patient. Returns the classical LoA of each bootstrap replicate.
bootstrap_observations <- function(d, n_boot, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  n <- length(d)

  # Column h holds the n measurements drawn in bootstrap replicate h
  boot_d <- matrix(sample(d, n * n_boot, replace = TRUE), nrow = n)
  bias <- colMeans(boot_d)
  sd_d <- sqrt((colSums(boot_d^2) - n * bias^2) / (n - 1))

  tibble(lower_loa = bias - z * sd_d,
         upper_loa = bias + z * sd_d)
}
