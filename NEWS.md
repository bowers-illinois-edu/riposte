# riposte 0.0.0.9005

* Polynomial rank scores as representations: `riposte_poly_scores(y, zeta)`
  gives the unit at within-block rank k of n_b the score
  `(k / (n_b + 1))^(zeta - 1)`, the score of Kim, Su, Bowers, and Li, and
  `riposte_poly_reps(zeta)` returns one representation per value of `zeta`.
  The default, `zeta = c(2, 7, 12, 17, 22)`, is the set in Bowers and Burton's
  rank-score tables. Unlike the Stephenson score, the polynomial score stays
  in [0, 1) in every block, and `zeta` need not be a whole number.
* The main vignette has a section on polynomial rank scores: how to supply
  your own `zeta`, how the largest `zeta` relates to the number of units per
  block whose gains the test should detect, and how to space the values so
  neighbouring scores are not near-copies of each other.
* The main vignette's statement about the default p-value now says what the
  exactness means: if the program changed no one's outcome, the chance of a
  p-value at or below 0.05 is at most 0.05, however many units or blocks there
  are.
* The covariance-adjustment paragraph of the main vignette no longer claims,
  without evidence in the package, that adjusting for covariates unrelated to
  the outcome loses almost no power. It now reports a simulation: with two covariates of pure noise and
  treatment raising the standard deviation from 1 to 1.35, the unadjusted test
  rejected in 85 of 100 experiments and the adjusted test in 83. That
  simulation is a new test in `tests/testthat/test-power.R`.

# riposte 0.0.0.9004

* `riposte_test()` accepts `engine = "asymptotic"` for unadjusted quadratic
  and Cauchy tests. The quadratic uses a chi-square distribution with degrees
  of freedom equal to the numerical rank used in its covariance inverse. The
  Cauchy combination uses the individual chi-square p-values and a standard
  Cauchy upper tail. These are approximations, not finite-sample exact tests.
  The permutation engine remains the default; no new runtime dependency is
  required.
* A formula without a block, `Y ~ treatment`, specifies complete randomization
  with the observed total treated count fixed. `Y ~ treatment | block` and
  the `blocks` argument continue to specify randomization within blocks.
* Asymptotic mode uses no random draws, ignores `seed` and `nresample`, and
  labels printed p-values as approximate. Both chi-square tails are calculated
  directly on the log scale for the Cauchy combination to avoid losing very
  small probabilities. An exactly zero component statistic gives a p-value of
  one and a negative infinite Cauchy term; it is not replaced by a permutation
  mid-p value. This differs from the paper simulation library's boundary
  fallback.
* The `riposte_acat_term()` examples and main vignette show the paper's Hybrid
  test: a Cauchy combination of the six individual p-values and the quadratic
  p-value, all at equal weight, using the existing functions.
* The asymptotic option requires an explicit `statistic = "quadratic"` or
  `"cauchy"`. Max, screen, and covariance-adjusted approximations are not
  implemented. Per-representation diagnostics from `riposte_components()`
  continue to use permutation mid-p values.

# riposte 0.0.0.9003

* The sixth default representation is now Huber's psi, `riposte_huber()`: the
  departure from the block median in units of the block's MAD, capped at
  +/- 1.345, dividing by the mean absolute departure when the MAD is zero. This
  matches the set of six in Bowers and Burton, "A More Powerful Test for
  Randomized Experiments" (working paper), and reproduces that paper's
  computation. The element of `riposte_reps_default()` is named `huber`;
  `riposte_tanh()` stays exported but is no longer in the default set.
  Results from `riposte_test()` and `riposte_components()` with the default
  representations change accordingly.
* `riposte_reps_proposal()` now swaps raw `tanh(y)` in for Huber's psi, so it
  still returns six representations.

# riposte 0.0.0.9001

* `riposte_test()` gains an optional `engine` argument. The default
  `engine = "permute"` is unchanged (exact, brute-force re-randomization).
  `engine = "saddlepoint"` computes the unadjusted Cauchy combination with no
  re-randomization: each representation's two-sided permutation p-value comes from
  a saddlepoint approximation to the exact within-block permutation distribution
  (via the optional `fastperm` package), and the Liu-Xie analytic tail combines
  them. It is an approximation, supports only `statistic = "cauchy"` without
  `adjust` (the quadratic and max need a multivariate saddlepoint not yet built),
  and is reached through `requireNamespace("fastperm")`, so it adds no formal
  dependency.

# riposte 0.0.0.9000

First working version. Block- and cluster-randomized combined randomization tests,
covariance adjustment, and the screen, all tested and `R CMD check`-clean.

## Methods

* Outcome representations as an open, user-supplied set
  (`riposte_reps_default()`, `riposte_reps_proposal()`, the individual
  `riposte_mean_dist()`, `riposte_mean_rank_dist()`, `riposte_max_dist()`,
  `riposte_rank()`, `riposte_tanh()`). The mean L1 distance is computed by a
  numerically stable centred closed form that matches the definitional
  sum-of-distances at any location offset (the manytestsr C++ prefix-sum form
  adds cancellation error there).
* Closed-form Strasser-Weber permutation moments (`riposte_sw_moments()`),
  verified against `coin` and against full enumeration.
* Three combinations on the representation linear statistics, all referred to the
  randomization distribution: the quadratic (full permutation covariance), the
  Cauchy (marginal mid-p, no clamp, no probit), and the max.
* The quadratic-vs-Cauchy screen (`riposte_screen()`), shrinking the covariance
  toward its diagonal by a conditioning-driven intensity; its condition number is
  exactly ancillary to the assignment.
* `riposte_test()` and `riposte_components()` --- the user-facing test and the
  per-representation diagnostic. Default statistic is the screen.
* Covariance adjustment (`riposte_adjust()`, `riposte_ridge_learner()`,
  `riposte_lm_learner()`): controls-only learner, residualize, refit inside every
  re-randomization (exact for any learner). The fit-once shortcut is available and
  documented as inexact under overfitting.
* Cluster-randomized designs (`clusters =`): collapse to the cluster level,
  permute clusters within blocks, effective n is the number of clusters.
* Combined Stephenson rank test as representations
  (`riposte_stephenson_reps()`, `riposte_stephenson_scores()`).

## Validation

* Size held under the sharp null for every combination (quadratic 0.050, Cauchy
  0.053, max 0.060) and for the screen and the cluster test.
* Power directions confirmed: the quadratic detects a cancelling effect and a
  scale change a difference in means misses; the Cauchy does not detect the
  cancelling effect; combined Stephenson beats a single Wilcoxon on a sparse
  upper-tail effect; ignoring clustering is anti-conservative.
* Exactness checked by full enumeration of tiny permutation distributions.
* Two rounds of adversarial multi-agent audit (correctness + the new
  cluster/Stephenson code), each finding verified by running R; all confirmed
  findings fixed.
* Policy-evaluator vignette (a Broader-Impacts deliverable).

## Exactness of the adjusted combinations

* The covariance-adjusted quadratic, max, and screen are exactly level-valid under
  the sharp null, for any learner and any `nresample`, when the metric is
  estimated from the POOLED set of statistics (observed + draws). Earlier code
  estimated it from the draws alone, which was anti-conservative (size ~0.12 at
  small `nresample`); the one-line fix pools. Proof and simulation evidence in
  `dev/theory-adjusted-exactness.md`; verified by simulation and by full
  enumeration. This resolves the two open theory items from the first handoff.

## Known limitations / open work

* The screen's shrinkage rule is a sensible placeholder, not the settled Project 1
  result; it is a power question only (exact level holds for any shrinkage rule).
* Cluster variance / few-cluster SE for ATE (Neyman) estimation is deferred to
  `propertee`; the Fisher (sharp-null) cluster test here is exact regardless.
* `propertee` interoperation is not yet implemented.
