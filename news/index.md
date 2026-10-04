# Changelog

## riposte 0.0.0.9007

- The quadratic combination no longer depends on the units of the
  outcome (issue
  [\#1](https://github.com/bowers-illinois-edu/riposte/issues/1)). It
  inverted the covariance of the score sums with
  [`MASS::ginv()`](https://rdrr.io/pkg/MASS/man/ginv.html), whose cutoff
  of sqrt(eps) times the largest eigenvalue dropped real directions when
  the rank sum’s variance was ~10^7 times the raw sum’s, as on survey
  scales; rescaling Y could then change Q and its p-value. It now
  divides each score sum by its standard deviation and inverts the
  correlation matrix (`riposte_std_pinv()`), which gives the same Q in
  exact arithmetic, and it reports as `df` the number of directions that
  inverse keeps (it reported [`qr()`](https://rdrr.io/r/base/qr.html)’s
  rank, which could exceed them). The asymptotic quadratic uses the same
  inverse. On an outcome with K \<= 7 values the quadratic now equals
  Pearson’s chi-square for the 2 x K table times (n - 1)/n on K - 1 df.
  Results change wherever `ginv()` had dropped a direction and are
  unchanged elsewhere.
- [`riposte_score_matrix()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_score_matrix.md)
  drops a representation as constant within blocks when its centred
  scores are within rounding error of zero relative to the size of its
  uncentred scores, instead of below the absolute number 1e-12, which
  dropped the raw score when Y was recorded in very small units.
- Two other fixed cutoffs of 1e-12 on a variance are gone for the same
  reason. Covariance adjustment drops a representation only when its
  statistics are constant relative to their own size, so it no longer
  drops the raw, distance, and max distance sums when Y is in small
  units. The check that refuses a design with no block holding both arms
  now looks for variances that are exactly zero, which is what such a
  design produces.

## riposte 0.0.0.9006

- [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
  and
  [`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
  take `alternative = "greater"` or `"less"` for a one-sided test of the
  sharp null of no effect. A one-sided test counts only treated-score
  sums above (or below) their re-randomization mean. It is available for
  the max combination, which becomes the default when `alternative` is
  one-sided, and for the Cauchy combination; the quadratic and the
  screen have no one-sided form. The permutation, asymptotic (Cauchy,
  from one-sided normal tails), and saddlepoint (Cauchy, from
  `fastperm`’s one-sided saddlepoint) engines and covariance adjustment
  all support it.
  [`riposte_midp()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_midp.md)
  takes the same argument.

- [`riposte_poly_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_scores.md)
  and
  [`riposte_poly_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_reps.md)
  take `tail = "lower"`, scores that weight the bottom of each block for
  a treatment that lowers the outcomes of a few units. The lower-tail
  score of the unit at rank k of n_b is
  `-((n_b + 1 - k) / (n_b + 1))^(zeta - 1)`; the minus sign makes it
  rise with the outcome, so `alternative = "less"` means lower treated
  outcomes for the raw outcome, the rank, and both tails.

- With rank-based representations, the default engine, no covariance
  adjustment, and no tied outcomes, the one-sided test is also valid for
  the hypothesis that the treatment lowered no unit’s outcome (Caughey,
  Dafoe, Li, and Miratrix 2023). A new test checks the inequality behind
  this on a full enumeration, and another checks that the one-sided max
  computes the same statistic as
  [`CMRSS::pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.html).

- Permutation p-values and mid-p values now count statistics within a
  relative tolerance of about 1.5e-8 as tied. Rank scores make many
  re-randomizations tie exactly, and floating-point summation separated
  such ties in their last bits, so whether a tied re-randomization
  counted as at least as extreme depended on rounding. Two exactness
  tests had matched that rounding: their designs (blocks of 4 with 2
  treated, one block of 4 clusters with 2 treated) give the same
  quadratic form at every assignment, so the exact p-value is 1. Those
  tests now use larger designs.

- A new vignette,
  [`vignette("harm")`](https://bowers-illinois-edu.github.io/riposte/articles/harm.md),
  looks for a few people harmed by a program, compares the one-sided
  test with the difference in means in a simulation, and uses `CMRSS`
  (now in Suggests) to ask how many were harmed and by how much. \#
  riposte 0.0.0.9005

- Polynomial rank scores as representations:
  `riposte_poly_scores(y, zeta)` gives the unit at within-block rank k
  of n_b the score `(k / (n_b + 1))^(zeta - 1)`, the score of Kim, Su,
  Bowers, and Li, and `riposte_poly_reps(zeta)` returns one
  representation per value of `zeta`. The default,
  `zeta = c(2, 7, 12, 17, 22)`, is the set in Bowers and Burton’s
  rank-score tables. Unlike the Stephenson score, the polynomial score
  stays in \[0, 1) in every block, and `zeta` need not be a whole
  number.

- The main vignette has a section on polynomial rank scores: how to
  supply your own `zeta`, how the largest `zeta` relates to the number
  of units per block whose gains the test should detect, and how to
  space the values so neighbouring scores are not near-copies of each
  other.

- The main vignette’s statement about the default p-value now says what
  the exactness means: if the program changed no one’s outcome, the
  chance of a p-value at or below 0.05 is at most 0.05, however many
  units or blocks there are.

- The covariance-adjustment paragraph of the main vignette no longer
  claims, without evidence in the package, that adjusting for covariates
  unrelated to the outcome loses almost no power. It now reports a
  simulation: with two covariates of pure noise and treatment raising
  the standard deviation from 1 to 1.35, the unadjusted test rejected in
  85 of 100 experiments and the adjusted test in 83. That simulation is
  a new test in `tests/testthat/test-power.R`.

## riposte 0.0.0.9004

- [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
  accepts `engine = "asymptotic"` for unadjusted quadratic and Cauchy
  tests. The quadratic uses a chi-square distribution with degrees of
  freedom equal to the numerical rank used in its covariance inverse.
  The Cauchy combination uses the individual chi-square p-values and a
  standard Cauchy upper tail. These are approximations, not
  finite-sample exact tests. The permutation engine remains the default;
  no new runtime dependency is required.
- A formula without a block, `Y ~ treatment`, specifies complete
  randomization with the observed total treated count fixed.
  `Y ~ treatment | block` and the `blocks` argument continue to specify
  randomization within blocks.
- Asymptotic mode uses no random draws, ignores `seed` and `nresample`,
  and labels printed p-values as approximate. Both chi-square tails are
  calculated directly on the log scale for the Cauchy combination to
  avoid losing very small probabilities. An exactly zero component
  statistic gives a p-value of one and a negative infinite Cauchy term;
  it is not replaced by a permutation mid-p value. This differs from the
  paper simulation library’s boundary fallback.
- The
  [`riposte_acat_term()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_acat_term.md)
  examples and main vignette show the paper’s Hybrid test: a Cauchy
  combination of the six individual p-values and the quadratic p-value,
  all at equal weight, using the existing functions.
- The asymptotic option requires an explicit `statistic = "quadratic"`
  or `"cauchy"`. Max, screen, and covariance-adjusted approximations are
  not implemented. Per-representation diagnostics from
  [`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
  continue to use permutation mid-p values.

## riposte 0.0.0.9003

- The sixth default representation is now Huber’s psi,
  [`riposte_huber()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_huber.md):
  the departure from the block median in units of the block’s MAD,
  capped at +/- 1.345, dividing by the mean absolute departure when the
  MAD is zero. This matches the set of six in Bowers and Burton, “A More
  Powerful Test for Randomized Experiments” (working paper), and
  reproduces that paper’s computation. The element of
  [`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md)
  is named `huber`;
  [`riposte_tanh()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_tanh.md)
  stays exported but is no longer in the default set. Results from
  [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
  and
  [`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
  with the default representations change accordingly.
- [`riposte_reps_proposal()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_proposal.md)
  now swaps raw `tanh(y)` in for Huber’s psi, so it still returns six
  representations.

## riposte 0.0.0.9001

- [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
  gains an optional `engine` argument. The default `engine = "permute"`
  is unchanged (exact, brute-force re-randomization).
  `engine = "saddlepoint"` computes the unadjusted Cauchy combination
  with no re-randomization: each representation’s two-sided permutation
  p-value comes from a saddlepoint approximation to the exact
  within-block permutation distribution (via the optional `fastperm`
  package), and the Liu-Xie analytic tail combines them. It is an
  approximation, supports only `statistic = "cauchy"` without `adjust`
  (the quadratic and max need a multivariate saddlepoint not yet built),
  and is reached through
  [`requireNamespace("fastperm")`](https://github.com/bowers-illinois-edu/fastperm),
  so it adds no formal dependency.

## riposte 0.0.0.9000

First working version. Block- and cluster-randomized combined
randomization tests, covariance adjustment, and the screen, all tested
and `R CMD check`-clean.

### Methods

- Outcome representations as an open, user-supplied set
  ([`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md),
  [`riposte_reps_proposal()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_proposal.md),
  the individual
  [`riposte_mean_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_mean_dist.md),
  [`riposte_mean_rank_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_mean_rank_dist.md),
  [`riposte_max_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_max_dist.md),
  [`riposte_rank()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_rank.md),
  [`riposte_tanh()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_tanh.md)).
  The mean L1 distance is computed by a numerically stable centred
  closed form that matches the definitional sum-of-distances at any
  location offset (the manytestsr C++ prefix-sum form adds cancellation
  error there).
- Closed-form Strasser-Weber permutation moments
  ([`riposte_sw_moments()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_sw_moments.md)),
  verified against `coin` and against full enumeration.
- Three combinations on the representation linear statistics, all
  referred to the randomization distribution: the quadratic (full
  permutation covariance), the Cauchy (marginal mid-p, no clamp, no
  probit), and the max.
- The quadratic-vs-Cauchy screen
  ([`riposte_screen()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_screen.md)),
  shrinking the covariance toward its diagonal by a conditioning-driven
  intensity; its condition number is exactly ancillary to the
  assignment.
- [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
  and
  [`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
  — the user-facing test and the per-representation diagnostic. Default
  statistic is the screen.
- Covariance adjustment
  ([`riposte_adjust()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_adjust.md),
  [`riposte_ridge_learner()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_ridge_learner.md),
  [`riposte_lm_learner()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_lm_learner.md)):
  controls-only learner, residualize, refit inside every
  re-randomization (exact for any learner). The fit-once shortcut is
  available and documented as inexact under overfitting.
- Cluster-randomized designs (`clusters =`): collapse to the cluster
  level, permute clusters within blocks, effective n is the number of
  clusters.
- Combined Stephenson rank test as representations
  ([`riposte_stephenson_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_reps.md),
  [`riposte_stephenson_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_scores.md)).

### Validation

- Size held under the sharp null for every combination (quadratic 0.050,
  Cauchy 0.053, max 0.060) and for the screen and the cluster test.
- Power directions confirmed: the quadratic detects a cancelling effect
  and a scale change a difference in means misses; the Cauchy does not
  detect the cancelling effect; combined Stephenson beats a single
  Wilcoxon on a sparse upper-tail effect; ignoring clustering is
  anti-conservative.
- Exactness checked by full enumeration of tiny permutation
  distributions.
- Two rounds of adversarial multi-agent audit (correctness + the new
  cluster/Stephenson code), each finding verified by running R; all
  confirmed findings fixed.
- Policy-evaluator vignette (a Broader-Impacts deliverable).

### Exactness of the adjusted combinations

- The covariance-adjusted quadratic, max, and screen are exactly
  level-valid under the sharp null, for any learner and any `nresample`,
  when the metric is estimated from the POOLED set of statistics
  (observed + draws). Earlier code estimated it from the draws alone,
  which was anti-conservative (size ~0.12 at small `nresample`); the
  one-line fix pools. Proof and simulation evidence in
  `dev/theory-adjusted-exactness.md`; verified by simulation and by full
  enumeration. This resolves the two open theory items from the first
  handoff.

### Known limitations / open work

- The screen’s shrinkage rule is a sensible placeholder, not the settled
  Project 1 result; it is a power question only (exact level holds for
  any shrinkage rule).
- Cluster variance / few-cluster SE for ATE (Neyman) estimation is
  deferred to `propertee`; the Fisher (sharp-null) cluster test here is
  exact regardless.
- `propertee` interoperation is not yet implemented.
