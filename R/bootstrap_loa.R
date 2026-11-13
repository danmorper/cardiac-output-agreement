library(tidyverse)

# Patient bootstrap for the limits of agreement.
# Resamples whole patients with replacement (keeping all their measurements),
# refits the model with fit_agreement() and returns the LoA of every refit.
# `response` is "difference" (additive model) or "log_ratio" (multiplicative model).
bootstrap_loa <- function(data, response, n_boot = 2000, seed = 1) {
  set.seed(seed)
  patients <- as.character(unique(data$patient))
  by_patient <- split(data, as.character(data$patient))

  map(seq_len(n_boot), function(b) {
    sampled <- sample(patients, replace = TRUE)

    # Each sampled copy gets a new id, so a patient drawn twice counts as two patients
    boot_data <- map2(by_patient[sampled], seq_along(sampled),
                      function(patient_data, k) mutate(patient_data, patient = paste0("B", k))) |>
      list_rbind()

    fit <- fit_agreement(boot_data, response = response)
    tibble(lower_loa = fit$lower_loa, upper_loa = fit$upper_loa)
  }) |>
    list_rbind()
}

# 95% percentile interval of a bootstrap distribution, as text "[a, b]"
percentile_ci <- function(x, digits = 2) {
  sprintf(paste0("[%.", digits, "f, %.", digits, "f]"),
          quantile(x, 0.025), quantile(x, 0.975))
}
