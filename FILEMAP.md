# File map: paper steps to code

Which file implements which step of the Estimation Details section, and which
file implements each bootstrap design. Intended as a checklist when verifying
the implementation against the paper.

Last updated after the bootstrap API was merged behind `bootstrap_ame()`.

---

## Estimation

The model is

```
Y[i,j,t] = a[i,t] + b[j,t] + x[i,j,t]'beta[t] + u[i,t]'v[j,t] + eps[i,j,t]
```

with each of the five parameter trajectories following a Gaussian random walk.
The estimator minimises the penalised MAP objective by block coordinate descent
at fixed penalties (inner loop), then refreshes the variance components by
empirical Bayes (outer loop).

| Paper step | File | Key function |
|---|---|---|
| Objective `Q` | `R/ame_objective.R` | `.ame_objective()` |
| **Step 1** initialization | `R/init_ame_trajectory.R` | `.init_ame_trajectory()` |
| **Step 2** update `beta` | `R/ame_blocks_beta.R` | `.update_beta_block()` |
| **Step 3** update `a` | `R/ame_blocks_additive.R` | `.update_a_block()` |
| **Step 4** update `b` | `R/ame_blocks_additive.R` | `.update_b_block()` |
| **Step 5** update `U` | `R/ame_blocks_latent.R` | `.update_U_block()` |
| **Step 6** update `V` | `R/ame_blocks_latent.R` | `.update_V_block()` |
| **Steps 7 + 8** identification | `R/ame_identify.R` | `.ame_identify()` |
| **Step 9** empirical-Bayes variances | `R/ame_eb_variances.R` | `.eb_update_variances()` |
| **Step 10** inner criterion | `R/ame_inner_bcd.R` | `.ame_inner_bcd()` |
| **Step 10** outer criterion | `R/ame_outer_eb.R` | `.ame_outer_eb()` |
| Numerical primitives | `R/tridiagonal_solvers.R` | block-tridiagonal solve, selected inversion, Procrustes, scale normalization |
| Public entry point | `R/fit_dynamic_ame.R` | `fit_dynamic_ame()` |

### Block structure of the trajectory updates

Each block update solves a block-tridiagonal system per independent unit. The
diagonal carries `information + lambda*[t=1] + gamma*(number of temporal
neighbours)`, i.e. `lambda + gamma` at `t = 1`, `2*gamma` in the interior,
`gamma` at `t = T`; the off-diagonals are `-gamma*I`.

| Block | Systems | Block size |
|---|---|---|
| `beta` | 1 (shared by every cell) | `P x P` |
| `a` | `N` (one per sender) | scalar |
| `b` | `M` (one per receiver) | scalar |
| `U` | `N` (one per sender) | `K x K` |
| `V` | `M` (one per receiver) | `K x K` |

**The right-hand side never contains a neighbour term.** In a joint-trajectory
solve the neighbours are unknowns, so the temporal coupling lives entirely in
the off-diagonal blocks. Adding `+gamma*neighbour` to the RHS is the sequential
filter formulation used by the deprecated `als_factorize_joint_cov()`, where the
previous period is fixed data; mixing it in here would count the coupling twice.

---

## Three deliberate departures from the paper as written

Each is documented in the header comment of the file that implements it.

### 1. Steps 7 and 8 are applied in the order 8 then 7

`R/ame_identify.R` normalises the scale first and rotates second, the reverse of
the paper's numbering.

Rotation preserves Frobenius norms, so applying it after the rescaling leaves
the Step-8 scale convention intact and both gauges hold at once. The reverse
does not work: rescaling changes the relative weight of `U` and `V` inside the
Procrustes objective, whose solution is driven by
`sum_t U_t'Uref + sum_t V_t'Vref`, which the rescaling turns into
`c*sum U_t'Uref + (1/c)*sum V_t'Vref`. Under rotate-then-rescale a second pass
returns a different rotation, so the reported latent coordinates keep turning
from one outer iteration to the next.

Both orders give identical fitted values and identical variance components —
only scale-then-rotate is idempotent (`R = I` to ~4e-16, `c = 1` exactly on a
second pass).

**The paper's Steps 7 and 8 should be swapped to match.**

### 2. Step 9 includes the posterior-variance term (a proper EM M-step)

`R/ame_eb_variances.R` adds trace corrections the paper's formulas omit:

```
E[ ||theta_1||^2 | Y ]           = ||theta_1_hat||^2 + tr(Sigma_11)
E[ ||theta_t - theta_{t-1}||^2 ] = ||d theta_hat||^2
                                   + tr(Sigma_tt) + tr(Sigma_{t-1,t-1})
                                   - 2 tr(Sigma_{t,t-1})
```

where `Sigma = sigma_eps^2 * H^{-1}` and `H` is the same block-tridiagonal
matrix the block updates solve. Only its diagonal and first sub-diagonal blocks
are needed; those come from the selected inversion in
`R/tridiagonal_solvers.R` at `O(T q^3)`, the same order as the solve.

Without these terms the plug-in update collapses: the point estimates are
already smoothed by the penalty, so their successive differences understate the
innovation, which raises the penalty, which smooths them further. Measured on a
20 x 14 x 8 panel, `tau_U2` fell from 9.5e-02 to the 1e-08 floor over 25 outer
iterations, `gamma_U` reached 4e7, and `cor(UV' hat, UV' true)` dropped to 0.24.
With the correction `tau_U2` settles at 7.3e-02 and the correlation is 0.97.

**The paper's Step 9 should add these terms.**

### 3. `tau_U2` and `tau_V2` are one parameter, not two

`R/ame_eb_variances.R` replaces both with their geometric mean
`kappa = sqrt(tau_U2 tau_V2)` before the penalties are derived, under
`latent_innovation = "shared"` (the default). `R/gibbs_dynamic_ame.R` does the
same for the sampler, drawing one `kappa` from the pooled sufficient statistics
of both sides rather than averaging two draws — the full conditional the model
implies, so the chain remains a correct sampler.

`U_t V_t'` is unchanged by `U_t -> c U_t`, `V_t -> V_t / c`, so the data
determine the product of the two innovation variances and not the split.
Estimated separately, the split is not measured but amplified: the update
derives `gamma_U = sigma_eps^2 / tau_U2` from the variance it has just
measured, so a smaller `tau_U2` smooths `U` further and lowers `tau_U2` again.

Measured, on 80 x 100 x 10 panels at `rho_x = 0.40`:

| generated | fitted `tau_U2 / tau_V2` | fitted product |
|---|---|---|
| ratio 4 (0.0098 / 0.00245) | 68.7, and still rising at pass 130 | within 20% |
| ratio 4, another panel | 49.8 | within 20% |
| ratio 1 (0.0049 / 0.0049) | 0.001 on one panel, 35 and 97 on others | within 20% |

The loop stopped in each case because the **fitted values** had settled: the
drift runs along a direction that leaves every `U_t V_t'` alone, so the
convergence criterion — a relative change in `Omega` and in `Yhat` — cannot see
it. Where the ratio lands is therefore a property of the stopping rule and not
of the data, which is why it must not be reported.

Both paths have to make the same change or the timing comparison in
simulations 3 and 4 is between two different models.

Two other repairs were tried first and are recorded because they look more
natural than the one that worked:

- **Per-period scale normalisation.** Choosing `c_t` from each period's
  magnitudes. Rejected on the argument that `c_t`'s own drift enters the
  increments as if it were temporal movement, and confirmed empirically: it
  reduced the fitted ratio from 39 to 19 and from 32 to 7.4 but never
  converged.
- **Anchoring the scale convention on the increments** instead of the
  magnitudes, so that `tau_U2 = tau_V2` holds after normalisation and the free
  coordinate moves to `sigma_V2 / sigma_U2`. This *breaks the estimator*: the
  factor does not return to one between passes, so the rescaling compounds —
  `c` sat at 3.3 for fifteen consecutive passes — and after forty passes the
  trajectories carry a factor near 1e7 and the solve goes singular. Reachable
  as `gauge = "innovation"`, kept only as a record.

The distinction that matters: the constraint acts on the **variance estimate**
and never touches the trajectories, so it has no mechanism by which to
compound. Imposing the same equality by rescaling does, and fails.

`latent_innovation = "separate"` restores the old behaviour, and
`fit$convergence$latent_ratio_trace` records what the unconstrained split would
have been at every pass whether or not it was used.

**The paper's Step 9 should estimate one latent innovation variance, and the
identification section should say why.**

---

## The outer loop can be accelerated

`R/ame_outer_eb.R`, `accelerate = "squarem"` (the default).

The outer loop is a fixed-point iteration on the ten penalties, and converges
linearly: on 80 x 100 x 10 it takes about 65 passes, while the inner loop it
drives settles in one or two. Nearly all of the cost of a fit is therefore in
how slowly the outer sequence approaches its limit.

A squared extrapolation (Varadhan and Roland 2008, package `SQUAREM`) over that
sequence reaches the same fixed point in fewer passes. Measured on
80 x 100 x 10:

```
plain      10.35 s   65 passes      one bootstrap refit  4.98 s
squarem     5.94 s   30 passes                           2.85 s
                                    largest beta difference  2.9e-04
```

Three things make it safe to have on by default:

- It works on the **log** penalties, so an extrapolated point cannot be
  negative however far it reaches.
- There is no objective to backtrack on — the penalised objective is not
  comparable across outer iterations, the penalties having changed — so after
  the extrapolation stops, the **plain iteration runs on until both original
  criteria hold**. What certifies the answer is the plain loop, in either mode.
- That certification gets the whole of `outer_max_iter`, not what the
  extrapolation left over. An earlier version shared one budget between the
  two; the extrapolation spends passes backtracking, and on two of eight real
  panels that left too few to finish, so the accelerated fit reported itself
  unconverged where the plain one had converged. An acceleration that converges
  worse than not accelerating is a defect, not a trade-off.

`accelerate = "none"` is the plain iteration and remains the definition of the
estimate. `tests/testthat/test-squarem.R` checks that the two reach the same
`Omega`, the same `beta`, and the same `U_t V_t'`.

---

## Bootstrap

`bootstrap_ame()` is the only exported bootstrap function; the two designs are
internal engines.

| Paper design | File | Function |
|---|---|---|
| Public entry point | `R/bootstrap_ame.R` | `bootstrap_ame()` |
| (a) conditional | `R/bootstrap_conditional.R` | `.bootstrap_conditional()` |
| (b) full-model | `R/bootstrap_full.R` | `.bootstrap_full()` |

```r
bt <- bootstrap_ame(fit, B = 1000, B_full = 500, n_cores = 4)
bt$conditional$summaries$beta     # estimate, SE, percentile CI
bt$conditional$latent$sign_prob   # Pr(U_t V_t' > 0) per dyad
bt$full$accuracy$Omega            # bias and RMSE against generating values
bt$full$latent_recovery           # Frobenius distance and correlation by period
```

### What each design reports, and why they differ

**Conditional.** Regenerates only `eps ~ N(0, sigma_eps2)` on top of the
systematic component. The dependence structure is preserved by construction
rather than by the resampling scheme, which is why the errors may be drawn
independently across cells. Reports standard errors and percentile intervals.

The systematic component it generates from is set by `dispersion`:

| `dispersion` | Generates from |
|---|---|
| `"posterior"` (default) | the fitted values with the dispersion the fitted `beta`, `a` and `b` are missing added back |
| `"none"` | the fitted values as they stand — the behaviour before the option existed |

A fitted trajectory is a posterior mean, and a posterior mean moves less than
the parameter it estimates. The empirical-Bayes update recovers `tau^2` only by
adding a posterior-variance trace to the squared increments of the fitted path,
so whatever that trace supplied is dispersion the path does not have. Measured
on the 80 x 100 x 5 simulation, the fitted `beta` path moves with innovation
variance 1.73e-04 against a `tau_beta^2` of 4.00e-04: **the fitted path is 66%
as wiggly as the truth.** Generating from it alone asks the estimator to track a
trajectory already as smooth as it wants trajectories to be, so the replicates
measure how the noise propagates but not the estimator's failure to follow a
moving truth. The reported standard errors came out **22% too small** (mean SE
over Monte Carlo SD = 0.78, uniformly across covariates and independent of
whether the fit converged).

`U` and `V` are never perturbed. `U_t V_t'` is unchanged by `(U_t A, V_t A^-T)`
for constant invertible `A`, and the identification step fixes only the
orthogonal and scalar parts of that freedom, so `tau_U^2` and `tau_V^2` can be
traded between the two factors without changing the fit — they describe the
chosen representative rather than the data. Injecting them as innovation noise
would create movement the data does not support. In the same simulation the
true `tau_U^2 = tau_V^2 = 4.9e-03` came back as 4.7e-03 and **2.14e-02**, a
factor of 4.4 apart for a symmetric truth, which is the gauge and not an error.

The correction is deliberately partial in a second way: the three blocks are
perturbed independently, while in the posterior they are correlated — they
compete for the same residual. Independent draws are more dispersed than joint
ones, and the ratio overshoots to about **1.6** as a result. A joint draw would
fix this, and `gibbs_dynamic_ame()` now provides one; see the note under the
sampler below.

**Full-model.** Redraws every trajectory from the zero-mean priors that the
fitted variance components parameterise, so each replicate has a *known*
generating value. Reports the difference between each estimate and the value
that generated it — mean error and RMSE — **not** the dispersion of the
estimates, which would be governed by how widely the generating values vary. An
estimator that recovered every replicate exactly would show the same dispersion.

### Conventions common to both designs

- Unobserved dyad-periods stay unobserved in every replicate.
- Covariates are held at their observed values; the model is conditional on them
  and does not model them.
- `K` is held fixed, so results do not reflect uncertainty about its selection.
- Every replicate re-estimates in full, including the variance-component update,
  so the penalties are re-derived within each replicate.
- Each replicate gets its own seed drawn up front, so results are identical
  regardless of core count or scheduling.

### Design-specific points

- The conditional design **warm-starts** each replicate at the original estimate
  (`Y*` is a small perturbation of `Y`, so this keeps replicates in one basin and
  makes the bootstrap measure sampling variability rather than optimiser
  variability). The full-model design **cannot** — its generating values differ —
  and runs the multistart in full, so it is several times slower per replicate.
- Exact percentile intervals for `U_t V_t'` require `store_latent_draws = TRUE`;
  otherwise a normal approximation is used and flagged in the returned object.
  The product's bootstrap distribution is skewed, so the two can differ
  materially.
- The full-model design reports `sigma_U2`, `sigma_V2`, `tau_U2`, `tau_V2` only
  through the gauge-invariant products `sigma_U2 * sigma_V2` and
  `tau_U2 * tau_V2`. The generated trajectories are drawn from the priors and do
  not satisfy the estimator's scale convention, whereas the estimates do, so the
  individual components would differ by an arbitrary gauge factor. Raw draws of
  all eleven are kept in `bt$full$draws$omega_full`.

### The interval the conditional design reports is a percentile interval

`.summarize_draws()` returns `quantile(beta*, 0.025)` and
`quantile(beta*, 0.975)`. That is centred on the bootstrap distribution, not
reflected through the estimate, and for a biased estimator the two differ in a
way that matters.

Under confounding the refits reproduce the original fit's bias, so `beta*` sits
off `beta_hat` in the same direction and by roughly the same amount as
`beta_hat` sits off `beta`. The percentile interval is then centred near
`beta_hat + bias` and **adds the bias a second time instead of removing it**.
The basic (reverse-percentile) interval

```
[ 2 * beta_hat - q_0.975(beta*),  2 * beta_hat - q_0.025(beta*) ]
```

inverts the error distribution the way the bootstrap principle intends, and can
be recomputed from what is already reported — no refitting needed.

Measured on the 80 x 100 x 5 simulation, coverage of a nominal 95% interval:

| | percentile | basic | offset of `beta*` from `beta_hat` |
|---|---|---|---|
| Static AME, `x3` (exogenous) | 0.988 | 0.989 | -0.01 |
| Static AME, `x1` (confounded with the latent structure) | 0.697 | 0.928 | 0.33 |
| Static AME, `x2` (confounded with the node effects) | **0.514** | **0.927** | 0.48 |

The offset is zero exactly where there is no confounding and largest where it
is strongest, which is what identifies the mechanism as the estimator's bias
rather than the interval arithmetic. **The package still reports percentile
limits**; the simulation scripts derive the basic interval from
`estimate`, `conf.low` and `conf.high` downstream. Switching the default is
open — see the list at the end of this file.

---

## Gibbs sampler

`R/gibbs_dynamic_ame.R` samples the posterior of the model
`fit_dynamic_ame()` maximises. It exists so the cost of the two inference paths
can be compared on equal terms — same likelihood, same priors, same data, same
machine, only the algorithm differing. The package's primary estimator remains
the penalised MAP fit.

| Purpose | Function | Exported |
|---|---|---|
| Draw from `N(H^-1 r, scale * H^-1)`, `H` block-tridiagonal | `.rtridiag_block()` | no |
| The same, vectorised over scalar systems | `.rtridiag_scalar_vec()` | no |
| One sweep over `beta, a, b, U, V` | `.gibbs_trajectories()` | no |
| The eleven variances | `.gibbs_omega()` | no |
| Entry point, on prepared arrays | `gibbs_dynamic_ame()` | **yes** |
| The same, from a `dynamic_ame` object | `gibbs_from_fit()` | **yes** |
| Posterior summaries, bulk and tail ESS | `summary.gibbs_dynamic_ame()` | **yes** |
| Rank-normalised split R-hat across chains | `gibbs_rhat()` | **yes** |

Every conditional is closed form, so there is no Metropolis step and nothing to
tune. **The sampler reuses the estimator's solvers rather than reimplementing
them**: `.solve_block_tridiagonal()`, `.block_tridiagonal_pivots()`,
`.solve_scalar_tridiagonal_vec()`, `.additive_diagonal()`,
`.latent_gram_rhs()`, `.cov_design()` and `.cov_term()` are called unchanged, so
no existing file was modified and the estimator carries no regression risk.

A sweep draws from the same block-tridiagonal system the block updates take the
mean of, at `O(T q^3)` — so **a sweep costs about one inner iteration of the
block coordinate descent**, which is the fact the timing comparison rests on.
The dispersion comes from the LDL' factorisation `H = L D' L'` whose pivots the
forward sweep already returns: with `L[t, t-1] = off * D'[t-1]^-1`, whitening by
`chol(D'[t])` and back-substituting through `L'` gives covariance `H^-1`.

### Two differences from the empirical-Bayes step

**No trace correction.** `.eb_update_variances()` adds a posterior-variance
trace to the squared increments because it plugs in a point estimate that is
itself smoothed. A Gibbs step conditions on a *drawn* trajectory, so no such
correction belongs there; adding it would count the same dispersion twice.

**The data arrive prepared.** The entry point takes `Y_list` and the covariate
lists, not an edge panel. Building those is shared work belonging to neither
inference path, and a timing comparison that charged it to one of them would be
measuring the wrong thing.

**One latent innovation variance.** `latent_innovation = "shared"` draws a
single `kappa` from the sufficient statistics of `U` and `V` pooled, not the
average of two draws — averaging two draws from separate conditionals is a
different distribution, and the chain would stop sampling the posterior it
claims to. This mirrors departure 3 in the estimator, and has to, or the timing
comparison is between two models.

### What to monitor

`U` and `V` are identified only up to a common orthogonal rotation and a
reciprocal rescaling, so their chains wander along that ridge indefinitely and
their R-hat never approaches one. That is a property of the parameter, not of
the sampler, and diagnosing on them would report a failure that is not there.
The reciprocal rescaling is fixed after every sweep — the scale gauge only;
`.scale_normalize_UV()` is called explicitly with `anchor = "pooled"` here
rather than taking the estimator's default, so that a change to that default
cannot move the sampler without anyone noticing. It did once.
Monitor `beta`, the variance components, or the products `U_t V_t'`. `beta` is
also what both inference paths report and identify, which makes it the right
basis for matching their precision.

`gibbs_rhat()` needs several chains started **far apart**: chains launched from
the same point agree immediately whether or not they have found the target. The
threshold in current use is 1.01 (Vehtari et al. 2021), not the older 1.1 — and
1.01 is what is consistent with an effective sample size in the hundreds, since
the residual between-chain discrepancy it allows is about 0.5% of a posterior
standard deviation against the 4.5% Monte Carlo error already accepted at
ESS 500, whereas at 1.1 the two would be the same size.

### It also makes an exact bootstrap possible

The conditional design's `dispersion = "posterior"` perturbs `beta`, `a` and `b`
independently, and overshoots because in the posterior they are correlated. A
draw from this sampler is a *joint* draw, so using one as the generating value
for each bootstrap replicate would fix that by construction. Not implemented and
not obviously worth implementing as a default — a user who is already running a
sampler can report the posterior interval directly — but valuable **once**, as a
way to measure what the independence approximation costs.

---

## Supporting files

| Purpose | File |
|---|---|
| Covariate term, design matrix, fitted values, residuals | `R/ame_helpers.R` |
| Panel matrix construction | `R/build_panel_matrices.R` |
| Covariate construction | `R/build_covariates.R` |
| Column validation | `R/checks.R` |
| `%||%` | `R/utils.R` |
| `print`, `summary`, `fitted`, `residuals`, and five accessors | `R/dynamic_ame_methods.R` |

---

## Deprecated, retained

The period-by-period estimator and everything built on it. These still run and
existing scripts are unaffected; only documentation marks them deprecated.

`als_factorize_joint_cov.R`, `fit_temporal_bipartite_als.R`,
`fit_temporal_bipartite_als_batch.R`, `fit_first_year_multistart.R`,
`svd_init.R`, `bootstrap.R`, `bootstrap_temporal_bipartite_als.R`

Supporting them: `accessors.R`, `extraction.R`, `summary_methods.R`, and
`post_estimation.R`. The last is also a compatibility layer for the new object —
`fit$results[[t]]` slices carry the legacy field names on purpose, so
`decompose_als_factorize_joint_cov()` and friends work on them directly.

`svd_init.R` is the only file that serves *nothing but* the deprecated path, and
would be the first to remove when the deprecated functions eventually go.

---

## Tests

| File | Covers |
|---|---|
| `test-tridiagonal_solvers.R` | solvers, selected inversion, Procrustes, scale normalization |
| `test-ame_objective.R` | `Q` against an element-by-element computation |
| `test-ame_blocks.R` | all five blocks against hand-built dense systems, exact-minimiser checks, monotonicity |
| `test-ame_eb_variances.R` | EM formulas, and a regression test for the variance collapse |
| `test-ame_identify.R` | invariance, idempotence, drift prevention |
| `test-fit_dynamic_ame.R` | inner loop, outer loop, public function, recovery |
| `test-dynamic_ame_methods.R` | methods and accessors |
| `test-bootstrap_ame.R` | the public wrapper |
| `test-bootstrap_conditional.R` | design (a), including the dispersion correction |
| `test-bootstrap_full.R` | design (b) |
| `test-gibbs_dynamic_ame.R` | the sampler |

Two tests in `test-gibbs_dynamic_ame.R` carry the sampler's correctness. The
first checks `.rtridiag_block()` against a densely constructed
`N(H^-1 r, scale * H^-1)` — mean within four Monte Carlo standard errors,
covariance within 6% over 20,000 draws. The second checks that the **posterior
mean of `beta` lands within half a posterior standard deviation of
`fit_dynamic_ame()`'s point estimate** on the same data. A sampler aimed at a
different target would pass the first and fail the second, which is what makes
it the load-bearing one: it is the evidence that the two inference paths are
inferring about the same model.

---

## Still to change in the paper

1. **Swap Steps 7 and 8** to match the implementation.
2. **Add the posterior-variance terms to Step 9.**
3. **State in both bootstrap designs** that unobserved dyad-periods remain
   unobserved.
4. The full-model design reports **seven gauge-free variance components plus two
   gauge-invariant products**, not all eleven.
5. Make the statement of what the latent factors are identified up to consistent
   between the two bootstrap sections: a common **orthogonal transformation**
   (sign changes are a special case) **and a reciprocal rescaling**.
6. Address the reviewer's remaining questions: how local minima are handled
   (warm start in the conditional design, full multistart in the full-model
   design, with the converged objective recorded as a diagnostic), and that `K`
   is held fixed.
7. **The conditional design generates from posterior means, which are smoother
   than the truth.** State the correction (`dispersion`), state that it is
   partial — `U` and `V` are left alone because `tau_U^2` and `tau_V^2` are
   gauge-dependent — and report that the independence of the three corrected
   blocks makes it overshoot. Measured: 0.78 before, about 1.6 after.
8. **Intervals should be basic, not percentile.** For a biased estimator the
   percentile interval adds the bias twice; Static AME coverage on `x2` goes
   from 0.514 to 0.927 under the reflected interval, while the exogenous
   covariate is unaffected. Either switch the package default or state in the
   paper which interval is reported.

---

## Open decisions in the package

- **Should `.summarize_draws()` report basic rather than percentile limits?**
  The evidence above says the percentile interval is wrong for this estimator,
  and the simulation scripts already correct it downstream. Changing the default
  alters results users may have saved, so it has not been done.
- **Should `dispersion` gain a joint-draw setting** built on
  `gibbs_dynamic_ame()`? See the task below; the answer depends on a
  measurement that has not been made.
- ~~**`tau_V2` came back 4.4x its generating value while `tau_U2` was
  accurate**~~ — **resolved.** The diagnosis was right and the consequence
  larger than it looked: on 80 x 100 panels the fitted ratio reached 68 and was
  still moving. The two are now one parameter; see departure 3.
- **Should the scale convention be anchored somewhere other than the pooled
  magnitudes?** `gauge = "first"` (the initial states) is available and is
  arguably the cleaner estimand — pooling lets a node set that genuinely moves
  more, and so spreads further, feed its own asymmetry back into the
  normalising constant. The default was left at `"pooled"` because the
  difference is numerically tiny once the innovation variances are tied, and
  because every result computed so far used it.
- **Should the sampler have been constrained at all?** The estimator's reason
  does not apply to it — a Gibbs step draws each variance from a full
  conditional rather than deriving a penalty from the value it just produced,
  so there is no feedback to run away. It was constrained for comparability,
  not because it was shown to drift, and nobody has checked whether it does.

---

## TODO: measure what the independence approximation costs

**Status: not done. Needed before the conditional bootstrap can be described
either way in the paper.**

### The problem

`dispersion = "posterior"` restores the missing dispersion to `beta`, `a` and
`b` **independently**. In the posterior those three are correlated — they
compete for the same residual, so a high draw of one goes with a low draw of
another. Independent draws are therefore more dispersed than joint ones, and
the correction overshoots:

| | mean SE / Monte Carlo SD |
|---|---|
| `dispersion = "none"` | 0.78 — too narrow |
| target | 1.00 |
| `dispersion = "posterior"` | **≈ 1.6 — too wide** |

Both numbers are wrong. Neither can be presented as the estimator's standard
error without saying which way it errs and by how much.

### The measurement

`gibbs_dynamic_ame()` produces draws from the *joint* posterior, which is
exactly what the approximation approximates. So:

1. Take a handful of the cached simulation fits.
2. Run the sampler on each to get `B` thinned joint draws
   `(beta, a, b, U, V)^(b)`.
3. For each draw, form `mu^(b)`, add `eps* ~ N(0, sigma_eps2)`, refit, collect
   `beta*` — the same conditional bootstrap, but generating from a joint draw
   instead of an independently perturbed one.
4. Compare the resulting `mean SE / Monte Carlo SD` against 0.78 and 1.6.

### How to read the answer

| Exact version lands at | Reading |
|---|---|
| ≈ 1.0 | The independence assumption is the whole error, and its cost is now quantified. Report the approximation with that number attached, or switch the default to the joint draw. |
| still ≫ 1 | The problem is **not** independence. Something else in the correction is wrong and has to be found before either version is used. |
| ≈ 0.78 | The dispersion correction itself is misconceived, not just its approximation. |

### Why it is worth doing even though users will not run it

Running an MCMC chain to repair a bootstrap is a strange thing to ask of a user
— someone already sampling can report the posterior interval and skip the
bootstrap entirely. This is not proposed as a default. It is a **one-off
validation**: it turns "we approximate the joint posterior by independent
perturbations" from an unexamined assumption into a statement with a measured
error attached, which is what the paper needs in order to say anything about
the interval at all.
