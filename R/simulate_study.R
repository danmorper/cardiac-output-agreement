# i=1,..., I patients, j=1,..., J_i repeated measurements per patient

# X_ij = true cardiac output for patient i at measurement j sampled from exp( [change between-patients] + [change within-patient] )
# change between-patients ~ N(log(co_median), co_sd_between) and change within-patient ~ N(0, co_sd_within)
# Y_ij = X_ij + [bias] + [change between-patients] + [change within-patient]
# D = Y_ij - X_ij = [bias] + [change between-patients] + [change within-patient] it is modeled as beta0 + u_i + e_ij, where u_i ~ N(0, sd_patient) and e_ij ~ N(0, sd_residual)

# Therefore I sample:
# - u_i ~ N(0, sd_patient) (change between-patients)
# - e_ij ~ N(0, sd_residual) (change within-patient)
# - X_ij ~ logN(log(co_median), co_sd_between) + N(0, co_sd_within) (true cardiac output)

# Y_ij is then computed as Y_ij = X_ij + bias + u_i + e_ij (additive model) or Y_ij = X_ij * exp(bias + u_i + e_ij) (multiplicative model)

simulate_study <- function(n_patients = 40,
                           n_per_patient = 8,
                           scale = c("additive", "multiplicative"),
                           bias = if (scale == "additive") 0.30 else log(1.05),
                           sd_patient = if (scale == "additive") 0.50 else 0.087,
                           sd_residual = if (scale == "additive") 0.70 else 0.123,
                           co_median = 5,
                           co_sd_between = 0.25,
                           co_sd_within = 0.15,
                           seed = NULL) {
  scale <- match.arg(scale)
  if (!is.null(seed)) set.seed(seed)

  # All patients have the same number of measurements
  if (length(n_per_patient) == 1) n_per_patient <- rep(n_per_patient, n_patients) # n_per_patient= 5 and n_patients = 6-> c(5, 5, 5, 5, 5, 5)

  # Basic checks
  if (any(n_per_patient < 1)) {
    stop("Every patient needs at least one measurement.")
  }
  if (sd_patient < 0 || sd_residual < 0) {
    stop("sd_patient and sd_residual must be non-negative.")
  }

  patient <- rep(seq_len(n_patients), times = n_per_patient)
  measurement <- sequence(n_per_patient)
  n <- length(patient)

  # Reference CO: log-normal around each patient's typical level.
  patient_log_co <- rnorm(n_patients, log(co_median), co_sd_between)
  reference <- exp(patient_log_co[patient] + rnorm(n, 0, co_sd_within))

  # Disagreement: random-intercept model.
  u <- rnorm(n_patients, 0, sd_patient)
  e <- rnorm(n, 0, sd_residual)
  disagreement <- bias + u[patient] + e

  test <- if (scale == "additive") reference + disagreement else reference * exp(disagreement)

  data <- data.frame(
    patient = factor(sprintf("P%03d", patient)),
    measurement = measurement,
    reference = reference,
    test = test,
    difference = test - reference,
    # The additive model can, very rarely, give a test value <= 0
    # (about 1 in 16,000 with the defaults). Its log ratio is then undefined.
    log_ratio = ifelse(test > 0, log(pmax(test, 1e-12) / reference), NA_real_),
    mean = (test + reference) / 2
  )

  attr(data, "truth") <- true_agreement(scale, bias, sd_patient, sd_residual)
  data
}

# Population bias, ICC and 95% limits of agreement implied by the parameters.
# For the multiplicative model the limits are ratios test / reference.
true_agreement <- function(scale, bias, sd_patient, sd_residual, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  sd_total <- sqrt(sd_patient^2 + sd_residual^2)
  lower <- bias - z * sd_total
  upper <- bias + z * sd_total

  out <- list(
    scale = scale,
    bias = bias,
    sd_patient = sd_patient,
    sd_residual = sd_residual,
    sd_total = sd_total,
    icc = sd_patient^2 / sd_total^2,
    lower_loa = lower,
    upper_loa = upper
  )
  if (scale == "multiplicative") {
    out$ratio_bias <- exp(bias)
    out$ratio_lower_loa <- exp(lower)
    out$ratio_upper_loa <- exp(upper)
  }
  out
}
