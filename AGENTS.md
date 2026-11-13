Write R code that is natural and easy to follow.

Use the tidyverse and ggplot2. If you use other libraries, explain why.

For tables in the Quarto documents use the gt package (`gt()`, `fmt_number()`, `cols_label()`), not printed tibbles or vectors.

Write in American English (analyze, modeled, summarize, color), also in code where both spellings exist (`summarize()`, `color =`).

Use simple and clear language, and avoid unnecessary complexity. Keep the code straightforward and easy to understand, but never at the cost of a mathematical or statistical mistake. This is a clinical study, at the level of an MSc in Mathematics.

## Project decisions

Already made. Keep them unless I say otherwise.

- **Names.** "Classical" analysis: Bland–Altman ignoring the patient (called naive in the literature). "Additive model": random intercept on $Y - X$. "Multiplicative model": random intercept on $\log(Y/X)$. In the thesis these were naive, M1 and M1log.
- **Confidence intervals.**
  - Bias: Wald interval. For the models, estimate ± $t_{0.975,\,I-1}$ · SE, with $I$ the number of patients (exact in a balanced design). For the classical analysis, $t_{0.975,\,n-1}$ with $n$ the number of measurements, as in Bland and Altman (1986).
  - Limits of agreement: patient (cluster) bootstrap with percentile intervals.
  - Written as `[a, b]`.
- **Caching.** Slow chunks use `cache: true` and `cache.extra` (see `analysis/agreement-analysis.qmd`).

### Monte Carlo simulation (`analysis/monte-carlo.qmd`)

- **Script.** `R/run_simulation.R` runs everything. I run it myself from the project folder (`Rscript R/run_simulation.R pilot`, then `Rscript R/run_simulation.R`). Do not run it unless I ask.
- **Output.** `results/study1.rds` and `results/study2.rds`, committed to git, one row per simulated dataset. The qmd only reads them: simulation chunks have `eval: false`, and all performance measures are computed in the qmd.
- **Sizes.** $B = 1{,}000$ datasets per scenario (Monte Carlo SE of a coverage near 95% is about 0.7%, Morris, White and Crowther 2019) and $H = 999$ bootstrap replicates per dataset.
- **Speed.** No `lmer()` in the loops. Every design is balanced, so the closed-form REML estimators of the one-way random-effects model are used (`R/fit_balanced.R`): $\hat\sigma^2_\varepsilon = MS_W$, $\hat\sigma^2_u = (MS_B - MS_W)/J$, and at the boundary ($MS_B < MS_W$) $\hat\sigma^2_u = 0$, $\hat\sigma^2_\varepsilon = SS_T/(n - 1)$. The bias SE is $\sqrt{MS_B/(IJ)}$ (SD of the patient means over $\sqrt I$), also at the boundary. The patient bootstrap resamples per-patient summaries. `monte-carlo.qmd` checks the closed form against `fit_agreement()`.
- **Study 1** (thesis): 4 ICC scenarios and 6 patients-versus-measurements scenarios. Classical analysis vs additive model. Three LoA CI methods: classical + observation bootstrap, classical + patient bootstrap, additive model + patient bootstrap.
- **Study 2** (thesis): multiplicative data, mild and marked scenarios, additive vs multiplicative model. The reference CO comes from `simulate_study()` (log-normal around 5 L/min), not the thesis' uniform 2–10 L/min.
- **Predictive coverage** is computed analytically from the true distribution, not by simulating new patients. In Study 2 it is a function of the reference CO.
- **Structure (decided 2026-11-08).** A short main thread plus an appendix. Main thread: intro, the headline figure "theory vs simulation" (predicted vs simulated coverage for every scenario) with the key findings, then three named ideas (1: design effect, 2: resample the independent unit, 3: constant limits vs an error that grows with CO), 6 figures and no tables. Full derivations, never simplified, go in `::: {.callout-note collapse="true" icon="false"}` boxes. Appendix A to E: Monte Carlo framework, closed-form REML, Study 1 in full, Study 2 in full (also the mild scenario), simulation code. All performance measures and predictions are computed in one folded `measures` chunk after the intro, because the headline figure needs them. The scenarios are read from `results/`, not from the script.
- **Theory first.** The results compare the simulation with the predictions of the ideas (predicted columns in tables, theory lines or crosses in figures). Study 1 builds on the independent $SS_W \sim \sigma_\varepsilon^2\chi^2_{I(J-1)}$, $SS_B \sim (J\sigma_u^2 + \sigma_\varepsilon^2)\chi^2_{I-1}$ and $\bar D$. Prediction 4 (observation bootstrap) is new, not in the thesis. Study 2 keeps the log-normal reference CO of `simulate_study()`: the limits the additive model converges to use log-normal moments, including $E(X_{ij}X_{ik}) = E(X)^2 e^{s_b^2}$ for two measurements of the same patient (decided 2026-11-08).
- **Reproducibility.** `parallel::mclapply()` over scenarios with `RNGkind("L'Ecuyer-CMRG")`, `set.seed()` and `mc.preschedule = FALSE`, so each scenario gets its own random-number stream whatever the number of cores.

### Joint agreement (`analysis/joint-agreement.qmd`)

- **Scale.** Additive model only, clinical parameters in L/min, as in the thesis application. On the log scale $\pi_J$, $R(J)$ and the tests would be the same, because $\exp()$ maps each coordinate one-to-one.
- **Data.** Two simulated devices with the same bias and total variance, so identical LoA: `simulate_study(scale = "additive", bias = 0.2, seed = 1)`, balanced 40 patients × 8 measurements, $V = \sigma_u^2 + \sigma_\varepsilon^2 = 0.40$. Device A: ICC 0.7 (patient-dominated). Device B: ICC 0.2 (residual-dominated). Only $J = 8$ is represented. Plots for $J = 2$ use measurements 1 and 2 of each patient.
- **Clinical parameters.** Thesis values $\beta_0^* = 0$, $\sigma_u^* = 0.65$, $\sigma_\varepsilon^* = 0.50$ L/min, 95% regions, $\alpha = 0.05$. Fixed before looking at the data.
- **Decision rules.** Primary: clinical coverage $\pi_J = P(\boldsymbol D_i \in A_{\mathrm{clin},J})$, acceptable if the one-sided 95% percentile patient-bootstrap lower bound is above 0.95 for every represented $J$. Secondary: containment, $R_{\max}$ with a one-sided 95% bootstrap upper bound below 1 (thesis). Both bounds are exploratory and not calibrated as tests, and the document says so. The exact binomial test (intersection–union over $J$) is shown with its limitation: it needs at least 59 patients with the same $J$.
- **Fitting.** Main fit with `fit_agreement()`. Bootstrap with `bootstrap_patients()` (closed-form REML, balanced design), 2,000 replicates.
- **Formulas.** Shared functions in `R/joint_agreement.R`. $\Sigma_J$ and $\Sigma_J^*$ share the eigenvectors $\boldsymbol 1_J/\sqrt J$ and its orthogonal complement, so $\pi_J$ is a 1-D integral and $R(J)$ has a closed form. Appendix E of the qmd compares them with Monte Carlo and numerical maximization.
- **Structure (decided 2026-11-09).** Same pattern as `monte-carlo.qmd`: a short main part plus an appendix. Main part: a short intro, the headline figure `fig-regions` (same LoA square, different ellipses, for $J = 2$) with three key findings about the ellipse, then five sections (same limits, different patients; the ellipse; the clinical region; counting patients is not enough; model-based criteria with the bootstrap), 6 figures and no tables. Full derivations, never simplified, go in `::: {.callout-note collapse="true" icon="false" title="Derivation"}` boxes. Appendix A to E: fitted parameters, the three regions for $J = 2$, the binomial test in full, the model-based criteria in full, checks. Data, the three regions, fits, criteria and the bootstrap are computed in one cached `measures` chunk after the intro, because the headline figure and the key findings need them. `fig-criteria-by-j` ($\pi_J$ and $R(J)$ against $J$) is in the model-based criteria section. Storyline (decided 2026-11-09): the ellipse is the main idea, the clinical criterion and the tests come after as its consequence.
- **No Monte Carlo** validation of the decision rules (decided 2026-11-08).

## Notes for the website version

The documents in `analysis/` will be turned into the project page https://danielmoreno.de/projects/cardiac-output-agreement.
Notes for that conversion are HTML comments that start with `WEB:`, for example `<!-- WEB: make this plot interactive -->`.
They do not show in the rendered document. Find them all with `grep -n "WEB:" analysis/*.qmd`.
When converting a document to the website, apply every `WEB:` note. Do not act on them while working on the analysis itself.

## Quarto: do not render the documents in `analysis/`

I usually keep `quarto preview` running on `analysis/*.qmd`, and it re-renders automatically every time I save. A second render of the same file (`quarto render`, another `quarto preview`, `Rscript -e "quarto::quarto_render(...)"`) writes to the same output and intermediate files (`.html`, `.knit.md`, `*_files/`) and breaks my preview.

- Never run `quarto render` or `quarto preview` on files in `analysis/`, even to check that your edit works.
- Never touch the documents' outputs (`.html`, `.knit.md`, `*_files/`, `*_cache/`).
- Do not start a preview yourself without asking me first.

### How to check your changes

1. Check whether a preview of the document you changed is running: `pgrep -fl "quarto preview"`.
2. **If it is, use it.** Note the modification time of the document's `.html`, save your edits, and wait until the preview has re-rendered (some documents take a few minutes, for example when a cache is invalidated).
   - **Newer `.html`: the render worked.** Check the result with `grep` on the `.html` (the text or numbers you expect), not by reading the whole page or taking screenshots.
   - **Same old `.html`: the render failed.** Quarto keeps the last good version, so the browser shows no error. If the preview runs in a terminal tab of the Claude app, read the error there. Otherwise, find the failing chunk by running the R code (next step).
3. **If no preview is running,** run only the R code (`source()` the scripts in `R/`, or run the chunks you changed with `Rscript`). If you really need a full render, copy the project to a scratchpad or temporary folder and render the copy there. Never in this folder.
4. Tell me what you ran and what you checked.

