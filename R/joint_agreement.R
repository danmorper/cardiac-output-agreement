library(tidyverse)

# Joint agreement for J repeated differences of one patient, used in analysis/joint-agreement.qmd.
#
# Under the additive model the patient difference vector D_i = (D_i1, ..., D_iJ) is
#   D_i ~ N_J(beta0 * 1_J, Sigma_J),   Sigma_J = sigma_u^2 * 1_J 1_J^T + sigma_e^2 * I_J.
# Sigma_J has two eigenvalues:
#   along 1_J / sqrt(J):                 lambda_J = J * sigma_u^2 + sigma_e^2
#   on the orthogonal complement (J - 1): sigma_e^2
# The clinical covariance Sigma*_J has the same eigenvectors, which is why the functions
# below reduce to one dimension.
#
# Clinical parameters are passed as a list: list(bias = , sd_patient = , sd_residual = ).

# Compound-symmetric covariance matrix Sigma_J
joint_covariance <- function(J, sd_patient, sd_residual) {
  sd_patient^2 * matrix(1, nrow = J, ncol = J) + sd_residual^2 * diag(J)
}

# TRUE for each patient (row of d) whose difference vector is inside the clinical
# region A_clin,J, the 95% ellipsoid of the clinical parameters
inside_clinical <- function(d, clinical, level = 0.95) {
  J <- ncol(d)
  distance <- mahalanobis(d, center = rep(clinical$bias, J),
                          cov = joint_covariance(J, clinical$sd_patient, clinical$sd_residual))
  distance <= qchisq(level, df = J)
}

# Points on the 95% ellipse for J = 2, for plotting.
# Its half-axes point along (1, 1) and (1, -1), with lengths radius * sqrt(eigenvalue).
ellipse_points <- function(bias, sd_patient, sd_residual, level = 0.95, n = 200) {
  radius <- sqrt(qchisq(level, df = 2))
  along <- radius * sqrt(2 * sd_patient^2 + sd_residual^2)  # along the diagonal (1, 1)
  across <- radius * sd_residual                            # across it, along (1, -1)
  angle <- seq(0, 2 * pi, length.out = n)
  tibble(d1 = bias + (along * cos(angle) + across * sin(angle)) / sqrt(2),
         d2 = bias + (along * cos(angle) - across * sin(angle)) / sqrt(2))
}

# Probability that all J differences are within half_width of the bias.
# Given the patient effect u, the J differences are independent N(bias + u, sigma_e^2),
# so the probability is a 1-D integral over u.
square_probability <- function(half_width, sd_patient, sd_residual, J = 2) {
  integrand <- function(u) {
    inside_one <- pnorm((half_width - u) / sd_residual) - pnorm((-half_width - u) / sd_residual)
    inside_one^J * dnorm(u, mean = 0, sd = sd_patient)
  }
  integrate(integrand, -Inf, Inf, rel.tol = 1e-8)$value
}

# Half-width of the square (cube for J > 2) that contains `level` joint probability.
# It lies between the LoA half-width (joint probability below `level`) and the
# Bonferroni half-width (joint probability at least `level`).
calibrated_half_width <- function(sd_patient, sd_residual, J = 2, level = 0.95) {
  alpha <- 1 - level
  sd_total <- sqrt(sd_patient^2 + sd_residual^2)
  uniroot(function(h) square_probability(h, sd_patient, sd_residual, J) - level,
          lower = qnorm(1 - alpha / 2) * sd_total,
          upper = qnorm(1 - alpha / (2 * J)) * sd_total,
          tol = 1e-10)$root
}

# Clinical coverage pi_J = P(D_i in A_clin,J) for a device with the given parameters.
#
# Split D_i - beta0* 1_J into a, its coordinate along 1_J / sqrt(J), and w, the part orthogonal to 1_J:
#   a ~ N(sqrt(J) * (beta0 - beta0*), lambda_J),   |w|^2 / sigma_e^2 ~ chi^2_{J-1},   independent.
# D_i is inside A_clin,J if a^2 / lambda*_J + |w|^2 / sigma*_e^2 <= chi^2_{J, level}, so
#   pi_J = integral over a of  density(a) * P(chi^2_{J-1} <= (sigma*_e / sigma_e)^2 * (cutoff - a^2 / lambda*_J)).
# For J = 1 there is no w: chi^2 with 0 degrees of freedom is 0, and pchisq() gives 1 inside the interval.
# One set of parameters at a time (use map or pmap for bootstrap replicates).
clinical_coverage <- function(J, bias, sd_patient, sd_residual, clinical, level = 0.95) {
  cutoff <- qchisq(level, df = J)
  lambda <- J * sd_patient^2 + sd_residual^2
  lambda_clinical <- J * clinical$sd_patient^2 + clinical$sd_residual^2
  shift <- sqrt(J) * (bias - clinical$bias)
  half_width <- sqrt(cutoff * lambda_clinical)  # largest |a| inside the clinical region

  integrand <- function(a) {
    room <- (cutoff - a^2 / lambda_clinical) * clinical$sd_residual^2 / sd_residual^2
    dnorm(a, mean = shift, sd = sqrt(lambda)) * pchisq(room, df = J - 1)
  }
  integrate(integrand, -half_width, half_width, rel.tol = 1e-8)$value
}

# Containment ratio R(J): the largest clinical distance (d - mu*)^T Sigma*_J^-1 (d - mu*) / cutoff
# over the device's 95% ellipsoid. R(J) <= 1 if and only if the device ellipsoid is inside A_clin,J.
#
# The distance is convex in d, so the maximum is on the boundary of the device ellipsoid.
# A boundary point has coordinate a = sqrt(cutoff * lambda_J) * t along 1_J / sqrt(J) and an
# orthogonal part with |w|^2 = cutoff * sigma_e^2 * (1 - t^2), for t in [-1, 1]. The distance
# depends on the point only through t, and it is a quadratic in t. So its maximum on [-1, 1]
# is at t = -1, t = 1 or at the vertex of the parabola.
# For J = 1 there is no orthogonal part and the boundary is just t = -1 and t = 1.
# Vectorized over bias, sd_patient and sd_residual (one value per bootstrap replicate).
containment_ratio <- function(J, bias, sd_patient, sd_residual, clinical, level = 0.95) {
  cutoff <- qchisq(level, df = J)
  lambda <- J * sd_patient^2 + sd_residual^2
  lambda_clinical <- J * clinical$sd_patient^2 + clinical$sd_residual^2
  shift <- sqrt(J) * (bias - clinical$bias)
  orthogonal <- if (J == 1) 0 else cutoff * sd_residual^2  # |w|^2 at t = 0

  clinical_distance <- function(t) {
    (sqrt(cutoff * lambda) * t + shift)^2 / lambda_clinical +
      orthogonal * (1 - t^2) / clinical$sd_residual^2
  }

  # Vertex of the parabola q2 * t^2 + q1 * t + q0, moved into [-1, 1]
  q2 <- cutoff * lambda / lambda_clinical - orthogonal / clinical$sd_residual^2
  q1 <- 2 * sqrt(cutoff * lambda) * shift / lambda_clinical
  vertex <- pmin(pmax(-q1 / (2 * q2), -1), 1)

  # na.rm: the vertex is undefined (0 / 0) when q2 = q1 = 0, then the distance is constant
  pmax(clinical_distance(-1), clinical_distance(1), clinical_distance(vertex), na.rm = TRUE) / cutoff
}
