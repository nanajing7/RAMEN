# RAMEN: Bipartite Dynamic Additive and Multiplicative Effects Network Model in R

## Acknowledgements

I thank Prof. Jared Edgerton and Prof. Ryan Kennedy for their invaluable guidance and support in the development of this package.

---

**RAMEN** fits dynamic additive and multiplicative effects (AME) models to bipartite panel data, where relationships between two sets of nodes evolve over time.

The model is

```
Y[i,j,t] = a[i,t] + b[j,t] + x[i,j,t]'beta[t] + u[i,t]'v[j,t] + eps[i,j,t]
```

with sender effects `a`, receiver effects `b`, regression coefficients `beta`, and latent factors `u`, `v`. Each of the five parameter trajectories follows a Gaussian random walk over time.

Typical applications include country–industry trade, legislator–bill co-sponsorship, and other bipartite relational panels.

---

## Estimation

`fit_dynamic_ame()` estimates all five trajectories **jointly across the whole panel**, rather than period by period.

Holding the other blocks fixed, each block's first-order conditions form a block-tridiagonal linear system — the random-walk prior couples each period only to its immediate neighbours — which is solved directly. Every such update is the exact minimiser of the objective over that block, so the objective is non-increasing across a sweep. The five blocks are swept until it stabilises.

The eleven variance components are then **estimated by empirical Bayes** rather than supplied as tuning constants, and the ten penalties are derived from them. The variance update is an EM M-step: it adds the posterior-variance term, without which the innovation variances collapse toward zero and take the fit with them.

Between the two loops the latent factors are identified — the scale gauge first, then a rotation onto a fixed reference — so that the variance components are well defined and the reported coordinates are reproducible.

```r
fit <- fit_dynamic_ame(
  edge_panel = my_panel,
  row_cov_df = my_row_covariates,
  row_covar_names = "gdp",
  K = 2
)

fit                       # dimensions, convergence, variance components
summary(fit)              # components alongside the penalties they imply
coef_trajectory(fit)      # beta_t, a P x T matrix
latent_trajectory(fit)    # U_t V_t' for each period
variance_components(fit)  # Omega, lambda, gamma
fitted(fit); residuals(fit); decompose_fit(fit)
```

Missing dyads are supported throughout: unobserved cells, and cells whose covariates are incomplete, are excluded from every sum rather than imputed.

---

## Inference

`bootstrap_ame()` runs two parametric bootstrap designs, which answer different questions. Both run by default and are returned together.

```r
bt <- bootstrap_ame(fit, B = 1000, B_full = 500, n_cores = 4)
bt                     # both designs side by side
```

**Conditional design — how precise are these estimates?** Holds the fitted systematic component fixed and regenerates only the observation errors, so every replicate carries exactly the dependence structure the model estimated. Gives standard errors and percentile confidence intervals for the parameters of the observed network.

```r
bt$conditional$summaries$beta    # estimate, std.error, conf.low, conf.high
bt$conditional$latent$sign_prob  # Pr(U_t V_t'[i,j] > 0) across replicates
```

**Full-model design — how well does the estimator work?** Redraws every state trajectory from the priors the fitted variance components parameterise and generates a complete new panel. Because the generating values are then known, the estimator is judged on the difference between each estimate and the value that produced it — bias and RMSE — rather than on the spread of its own output, which would reflect the variation of the generating values instead.

```r
bt$full$accuracy$Omega     # mean error and RMSE against the generating values
bt$full$latent_recovery    # Frobenius distance and correlation, by period
```

Either can be run alone with `design = "conditional"` or `design = "full"`.

Replications are independent and parallelise cleanly. Each is given its own seed up front, so results do not depend on the number of workers or on scheduling.

### What is identified

`U_t V_t'` is identified; the individual factors are not, beyond the gauge the fit fixes. Interpret rotation-invariant quantities — the fitted multiplicative component, the fitted values, the coefficient trajectories. For the same reason `sigma_U2`, `sigma_V2`, `tau_U2`, and `tau_V2` are reported in the full-model bootstrap only through their gauge-invariant products.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("nanajing7/RAMEN")
```

## Demo

`demo/demo_dynamic_ame.R` simulates a small panel, fits it, and shows the accessors and both bootstrap designs. It runs in well under a minute.

## Paper replication

The `paper_example_simulation/` folder contains the simulation scripts used to reproduce the results reported in the accompanying paper.

## Deprecated

The earlier period-by-period alternating least squares functions — `fit_temporal_bipartite_als()`, `fit_temporal_bipartite_als_batch()`, `als_factorize_joint_cov()`, `fit_first_period_multistart()`, and the two bootstrap functions built on them — still run, and existing scripts are unaffected. They fit each period in a forward sweep with the previous period held fixed, and take the penalties as user-supplied constants. Use `fit_dynamic_ame()` for new work.
