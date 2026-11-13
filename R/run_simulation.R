# Monte Carlo studies for analysis/monte-carlo.qmd.
#
# Run from the project folder (cardiac-output-agreement/):
#   Rscript R/run_simulation.R pilot   # B = 20: prints the run time and an estimate for B = 1,000. Saves nothing.
#   Rscript R/run_simulation.R         # B = 1,000: saves results/study1.rds and results/study2.rds
#
# Each .rds has one row per simulated dataset: the scenario, the estimates of every
# method and the 95% bootstrap CIs of the limits of agreement. The performance
# measures (coverage, bias, ...) are computed in analysis/monte-carlo.qmd.
#
# The lines "## ---- name ----" mark the parts shown in analysis/monte-carlo.qmd.

library(tidyverse)
library(parallel)  # mclapply() runs the scenarios in parallel, part of base R

source("R/simulate_study.R")
source("R/fit_balanced.R")

pilot <- "pilot" %in% commandArgs(trailingOnly = TRUE)
n_sim <- if (pilot) 20 else 1000  # B: simulated datasets per scenario
n_boot <- 999                     # H: bootstrap replicates per dataset
n_cores <- 8

## ---- study1-scenarios ----
# The bias and the total variance V of D are fixed. Only their split between patients
# (ICC) changes, so every scenario has the same true bias and limits of agreement.
V <- 0.74

icc_scenarios <- tibble(
  scenario = c("No clustering", "Residual-dominated", "Moderate clustering", "Patient-dominated"),
  block = "ICC",
  icc = c(0, 0.100, 0.338, 0.750),
  n_patients = 40,
  n_per_patient = 8
)

# Same 160 measurements, split between patients in different ways
allocation_scenarios <- tibble(
  block = "Allocation",
  icc = rep(c(0.338, 0.750), each = 3),
  n_patients = rep(c(20, 40, 80), times = 2),
  n_per_patient = rep(c(8, 4, 2), times = 2)
) |>
  mutate(scenario = paste0("ICC ", icc, ": ", n_patients, " patients x ", n_per_patient))

scenarios_1 <- bind_rows(icc_scenarios, allocation_scenarios) |>
  mutate(true_bias = 0.30,
         sd_patient = sqrt(icc * V),
         sd_residual = sqrt((1 - icc) * V))

## ---- study2-scenarios ----
# Multiplicative data: geometric mean ratio Y/X of 1.05 and log-scale ICC of 0.338.
# Only the total SD of log(Y/X) changes.
scenarios_2 <- tibble(
  scenario = c("Mild", "Marked"),
  sd_log_total = c(0.10, 0.20),
  icc = 0.338,
  n_patients = 40,
  n_per_patient = 8
) |>
  mutate(true_bias = log(1.05),
         sd_patient = sqrt(icc) * sd_log_total,
         sd_residual = sqrt(1 - icc) * sd_log_total)

## ---- helpers ----
# 95% percentile CIs of the lower and upper LoA from their bootstrap replicates,
# as one row with the column names prefixed by the method
loa_ci <- function(boot_lower, boot_upper, method) {
  tibble(lower_loa_ci_low  = quantile(boot_lower, 0.025, names = FALSE),
         lower_loa_ci_high = quantile(boot_lower, 0.975, names = FALSE),
         upper_loa_ci_low  = quantile(boot_upper, 0.025, names = FALSE),
         upper_loa_ci_high = quantile(boot_upper, 0.975, names = FALSE)) |>
    rename_with(\(name) paste0(method, "_", name))
}

# Estimates of the random-intercept model, with the column names prefixed by the model
model_columns <- function(fit, model) {
  fit |>
    select(bias, se, sd_patient, sd_residual, boundary, lower_loa, upper_loa) |>
    rename_with(\(name) paste0(model, "_", name))
}

## ---- study1-analysis ----
# Study 1: simulate one dataset and analyze it with the classical analysis and the additive model
analyze_study1 <- function(scenario) {
  data <- simulate_study(n_patients = scenario$n_patients,
                         n_per_patient = scenario$n_per_patient,
                         scale = "additive",
                         bias = scenario$true_bias,
                         sd_patient = scenario$sd_patient,
                         sd_residual = scenario$sd_residual)

  # One row per patient, one column per measurement (the data are sorted by patient)
  d <- matrix(data$difference, nrow = scenario$n_patients, byrow = TRUE)

  fit <- fit_balanced(d)
  boot_obs <- bootstrap_observations(d, n_boot)  # resamples single measurements
  boot_pat <- bootstrap_patients(d, n_boot)      # resamples whole patients

  bind_cols(
    select(fit, starts_with("classical_")),
    model_columns(fit, "additive"),
    loa_ci(boot_obs$lower_loa, boot_obs$upper_loa, "classical_obs"),
    loa_ci(boot_pat$classical_lower_loa, boot_pat$classical_upper_loa, "classical_patient"),
    loa_ci(boot_pat$lower_loa, boot_pat$upper_loa, "additive_patient")
  )
}

## ---- study2-analysis ----
# Study 2: simulate one multiplicative dataset and fit the additive model (to Y - X)
# and the multiplicative model (to log(Y / X)), each with its patient bootstrap
analyze_study2 <- function(scenario) {
  data <- simulate_study(n_patients = scenario$n_patients,
                         n_per_patient = scenario$n_per_patient,
                         scale = "multiplicative",
                         bias = scenario$true_bias,
                         sd_patient = scenario$sd_patient,
                         sd_residual = scenario$sd_residual)

  d_raw <- matrix(data$difference, nrow = scenario$n_patients, byrow = TRUE)
  d_log <- matrix(data$log_ratio, nrow = scenario$n_patients, byrow = TRUE)

  boot_add <- bootstrap_patients(d_raw, n_boot)
  boot_mult <- bootstrap_patients(d_log, n_boot)

  # The multiplicative estimates and CIs stay on the log scale (exponentiate for ratios)
  bind_cols(
    model_columns(fit_balanced(d_raw), "additive"),
    model_columns(fit_balanced(d_log), "multiplicative"),
    loa_ci(boot_add$lower_loa, boot_add$upper_loa, "additive"),
    loa_ci(boot_mult$lower_loa, boot_mult$upper_loa, "multiplicative")
  )
}

## ---- run ----
# All n_sim datasets of one scenario, one row per dataset
run_scenario <- function(scenario, analyze) {
  start <- Sys.time()
  results <- map(seq_len(n_sim), \(b) analyze(scenario)) |>
    list_rbind() |>
    mutate(replicate = seq_len(n_sim), n_boot = n_boot, .before = 1)

  seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  message(sprintf("  %-38s %6.1f s", scenario$scenario, seconds))
  cross_join(scenario, results)
}

# One parallel job per scenario. With mc.preschedule = FALSE, scenario k always gets
# the k-th random-number stream, whatever the number of cores.
run_study <- function(scenarios, analyze) {
  results <- mclapply(seq_len(nrow(scenarios)),
                      \(k) run_scenario(scenarios[k, ], analyze),
                      mc.cores = n_cores, mc.preschedule = FALSE)

  # mclapply() returns an error object instead of stopping when a job fails
  failed <- map_lgl(results, \(x) inherits(x, "try-error"))
  if (any(failed)) stop("A scenario failed: ", results[failed][[1]])
  list_rbind(results)
}

RNGkind("L'Ecuyer-CMRG")  # random-number generator with independent parallel streams
start <- Sys.time()

message(sprintf("Study 1: %d scenarios, B = %d, H = %d", nrow(scenarios_1), n_sim, n_boot))
set.seed(1)
study1 <- run_study(scenarios_1, analyze_study1)

message(sprintf("Study 2: %d scenarios, B = %d, H = %d", nrow(scenarios_2), n_sim, n_boot))
set.seed(2)
study2 <- run_study(scenarios_2, analyze_study2)

minutes <- as.numeric(difftime(Sys.time(), start, units = "mins"))

if (pilot) {
  message(sprintf("Pilot (B = %d) took %.1f min. Full run (B = 1,000): about %.0f min.",
                  n_sim, minutes, minutes * 1000 / n_sim))
} else {
  dir.create("results", showWarnings = FALSE)
  write_rds(study1, "results/study1.rds")
  write_rds(study2, "results/study2.rds")
  message(sprintf("Done in %.1f min. Saved results/study1.rds and results/study2.rds.", minutes))
}
