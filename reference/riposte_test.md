# Combined randomization-based test of the sharp null

Combined randomization-based test of the sharp null

## Usage

``` r
riposte_test(
  formula,
  data,
  blocks = NULL,
  clusters = NULL,
  statistic = c("screen", "quadratic", "cauchy", "max"),
  representations = riposte_reps_default(),
  nresample = 1999L,
  screen = riposte_screen_control(),
  adjust = NULL,
  learner = NULL,
  adjust_refit = TRUE,
  cluster_agg = mean,
  seed = NULL,
  engine = c("permute", "saddlepoint", "asymptotic"),
  alternative = c("two.sided", "greater", "less"),
  ...
)
```

## Arguments

- formula:

  a design formula `Y ~ treatment` for complete randomization, or
  `Y ~ treatment | block` for randomization within blocks (or use
  `blocks`).

- data:

  a data.frame/data.table.

- blocks:

  optional block specification: a vector of labels or a column name in
  `data`. Omit when the block is in the formula as `| block`. With
  neither specification, the design fixes the total treated count.

- clusters:

  optional cluster specification for a cluster-randomized design: a
  vector of cluster labels or a column name in `data`. When given, the
  design is collapsed to one row per cluster (outcome and covariates
  aggregated by `cluster_agg`), treatment is permuted over clusters
  within blocks, and the effective sample size is the number of
  clusters. Treatment must be constant within a cluster and clusters
  must nest within blocks.

- statistic:

  which combination to use: `"screen"` (default), `"quadratic"`,
  `"cauchy"`, or `"max"`. With a one-sided `alternative`, only `"max"`
  (the default then) and `"cauchy"` are available.

- representations:

  a named list of representation functions;
  [`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md)
  by default. Supply your own to combine a different set.

- nresample:

  number of re-randomizations for the Monte-Carlo combinations. Ignored
  for `engine = "asymptotic"`.

- screen:

  a
  [`riposte_screen_control()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_screen_control.md)
  list, used when `statistic = "screen"`.

- adjust:

  optional covariance adjustment: either a one-sided covariate formula
  (e.g. `~ x1 + x2`) or a
  [`riposte_adjust()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_adjust.md)
  specification. When a bare formula, it is wrapped with `learner`
  (ridge by default). `NULL` means no adjustment. The combination is
  then computed on the controls-only residuals, refit inside each
  re-randomization; see
  [`riposte_adjust()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_adjust.md)
  for what is and is not exact under adjustment.

- learner:

  optional controls-only learner `(X, Y, control) -> Yhat` for
  covariance adjustment;
  [`riposte_ridge_learner()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_ridge_learner.md)
  by default. Overrides the learner in an `adjust` specification if both
  are given.

- adjust_refit:

  refit the learner inside every re-randomization (`TRUE`, exact for any
  learner) or fit once and permute fixed residuals (`FALSE`, the
  shortcut). Used only when `adjust` is a bare formula.

- cluster_agg:

  how to aggregate the outcome (and covariates) within a cluster when
  `clusters` is given; the cluster mean by default.

- seed:

  optional integer seed (L'Ecuyer-CMRG) for reproducibility. Ignored for
  `engine = "asymptotic"`, which uses no random draws.

- engine:

  `"permute"` (default) refers every combination to brute-force
  re-randomization (exact). `"asymptotic"` uses a chi-square reference
  for `statistic = "quadratic"`, with degrees of freedom equal to the
  numerical rank of the permutation covariance. For
  `statistic = "cauchy"`, it computes each representation's two-sided
  p-value from its squared standardized statistic and a chi-square
  distribution on one degree of freedom, then compares the mean Cauchy
  transform with a standard Cauchy distribution. These are
  approximations, available without `adjust`, and do not require `coin`
  or `fastperm`. Set `statistic` explicitly to `"quadratic"` or
  `"cauchy"`; the max and screen approximations are not implemented.
  `"saddlepoint"` computes the combination with no re-randomization, via
  `fastperm`: the unadjusted Cauchy (each representation's marginal
  saddlepoint p-value combined by the Liu-Xie analytic tail), the
  unadjusted quadratic (the omnibus `Q` referred to fastperm's
  multivariate saddlepoint M2, feasible for a few representations), or
  `"screen"`, which uses the permutation `Sigma` to choose between them
  (and to regularize a rank-deficient metric). It is an approximation,
  supports `statistic` other than `"max"` without `adjust`, and requires
  the `fastperm` package (the quadratic needs a `fastperm` with
  `fastperm_spa_quadratic`).

- alternative:

  `"two.sided"` (default) counts evidence in either direction.
  `"greater"` counts only treated-score sums above their
  re-randomization mean, and `"less"` only sums below it. For scores
  that rise with the outcome — the raw outcome, the rank, and both tails
  of
  [`riposte_poly_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_reps.md)
  — `"less"` asks whether treated outcomes are lower, which is how to
  look for harm. The distance representations do not rise with the
  outcome, so a one-sided test of them asks instead whether treated
  units sit closer to (`"less"`) or farther from (`"greater"`) the rest
  of their block. Every alternative tests the same sharp null hypothesis
  of no effect; see Details.

- ...:

  reserved.

## Value

an object with the observed statistic, the p-value, the combination
used, the condition number, and (for `"screen"`) the choice made. With
`engine = "permute"`, the p-value is the share of the `nresample`
re-randomizations at least as extreme as observed, with the observed
assignment included. That observed-inclusion convention makes the test's
level exactly valid in finite samples for any `nresample` (the p-value
estimates, rather than equals, the full-enumeration value, so it varies
with the seed). The closed-form moments and the screen's condition
number are exact, not Monte-Carlo. With `engine = "asymptotic"`,
`nresample` is zero and the p-value uses the stated approximation, not a
permutation count. Cauchy results from this engine also include
`component_p`, the individual chi-square p-values.

## Details

The chi-square approximation requires a normal limit for the
standardized treated-score sums on their nondegenerate subspace. A large
number of units alone does not ensure this when individual scores
dominate their variance. For the Cauchy combination, dependent component
p-values do not in general yield an exactly standard Cauchy statistic:
the analytic upper tail is an approximation, not a finite-sample level
guarantee.

A one-sided test is still a test of the sharp null of no effect, and it
holds its level under that null like the two-sided test. With `"less"`,
the default engine, no covariance adjustment, representations that
depend on the outcome only through its within-block ranks, and no tied
outcomes, it also holds its level under the weaker hypothesis that the
treatment lowered no unit's outcome, though it may have raised some
(Caughey, Dafoe, Li, and Miratrix 2023, "Randomisation inference beyond
the sharp null", JRSS-B 85: 1471-1491). The ranks in each block are the
same numbers whatever the treatment did, so the re-randomization
distribution does not change, and raising treated outcomes can only
raise the observed treated-score sums, so the p-value can only grow. The
same holds for `"greater"` with "raised" and "lowered" exchanged.

The asymptotic Cauchy calculation retains both tails on the log scale to
avoid artificial zeros and ones. An exactly zero component statistic
gives an individual p-value of one and a Cauchy term of negative
infinity; the combined p-value is then one unless an opposing infinite
term makes the sum undefined. That case gives an error. This engine
never switches silently to permutations. It combines the individual
p-values only; it does not add the quadratic p-value as another term.
The example in
[`riposte_acat_term()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_acat_term.md)
shows how to include that p-value to calculate the paper's Hybrid test.
