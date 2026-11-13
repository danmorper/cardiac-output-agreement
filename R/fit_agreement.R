library(tidyverse)
library(lme4)  # linear mixed models with lmer(), not part of the tidyverse

# Fits the random-intercept model D_ij = beta0 + u_i + e_ij by REML and
# returns bias (with 95% Wald CI), variance components and 95% limits of agreement.
# `response` is "difference" for the additive model or "log_ratio" for the multiplicative model.
fit_agreement <- function(data, response = c("difference", "log_ratio"), level = 0.95) {
  response <- match.arg(response)
  z <- qnorm(1 - (1 - level) / 2)

  # The bias is estimated from I independent patients, so its CI uses a t quantile
  # with I - 1 degrees of freedom (exact in a balanced design, unless sd_patient = 0).
  n_patients <- n_distinct(data$patient)
  t_crit <- qt(1 - (1 - level) / 2, df = n_patients - 1)

  fit <- lmer(reformulate("1 + (1 | patient)", response = response),
              data = data, REML = TRUE)

  bias <- fixef(fit)[[1]]
  se_bias <- sqrt(vcov(fit)[1, 1])

  # VarCorr rows: first the patient random intercept, then the residual
  sds <- as.data.frame(VarCorr(fit))$sdcor
  sd_patient <- sds[1]
  sd_residual <- sds[2]
  sd_total <- sqrt(sd_patient^2 + sd_residual^2)

  out <- tibble(
    model = if (response == "difference") "additive" else "multiplicative",
    bias = bias,
    bias_lower = bias - t_crit * se_bias,
    bias_upper = bias + t_crit * se_bias,
    sd_patient = sd_patient,
    sd_residual = sd_residual,
    icc = sd_patient^2 / sd_total^2,
    lower_loa = bias - z * sd_total,
    upper_loa = bias + z * sd_total
  )

  # On the log scale, exponentiate to get limits for the ratio Y / X
  if (response == "log_ratio") {
    out <- out |>
      mutate(ratio_bias = exp(bias),
             ratio_lower_loa = exp(lower_loa),
             ratio_upper_loa = exp(upper_loa))
  }
  out
}
