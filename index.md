# riposte

Tests that combine several representations of an outcome to detect
treatment effects a difference of means misses, in completely
randomized, block-randomized, and cluster-randomized experiments. Users
can choose direct re-randomization or, for the unadjusted quadratic and
Cauchy tests, an analytic approximation.

The name: Randomization Inference for Power via Outcome representationS,
Tests, and Estimation.

## Installation

riposte is not on CRAN. Install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("bowers-illinois-edu/riposte")
```

The optional saddlepoint engine also needs the `fastperm` package:

``` r

remotes::install_github("bowers-illinois-edu/fastperm")
```

## What it does

Given an outcome Y and a block- (or cluster-) randomized design, riposte
builds several within-block representations of Y — raw values, ranks,
mean pairwise distance, mean rank distance, max distance, and a bounded
transform — and combines the evidence each carries:

- a **quadratic** (energy/distance) omnibus that uses the full
  permutation covariance of the representations;
- a **Cauchy** combination that uses only each representation’s marginal
  calibration; and
- a **max** combination, the largest standardized representation; and
- a **min-p** combination (`statistic = "minp"`), the smallest of the
  representations’ mid-p values, which equals the max when every
  representation’s statistic has the same null distribution.

A **screen** chooses between the quadratic and the Cauchy from the
conditioning of the permutation covariance — an ancillary quantity (a
function of the outcomes under the sharp null and the design, not of the
realized assignment), so choosing from it does not disturb the exact
level.

By default, p-values use direct re-randomization under the specified
design. `engine = "asymptotic"` instead uses a chi-square reference for
the quadratic statistic. For the Cauchy combination, it computes
individual chi-square p-values and compares their mean Cauchy transform
with a standard Cauchy distribution. Both calculations are
approximations; they do not have the permutation test’s finite-sample
level guarantee.

The examples below use a simulated experiment: 20 blocks, each with six
clusters of two units, and three clusters treated in each block.
Treatment doubles the standard deviation of the outcome and leaves its
mean unchanged, an effect a difference of means cannot detect.

``` r

library(riposte)
set.seed(20260929)
d <- data.frame(block   = factor(rep(1:20, each = 12)),
                cl      = factor(paste(rep(1:20, each = 12),
                                       rep(rep(1:6, each = 2), times = 20))),
                treated = rep(rep(0:1, each = 6), times = 20),
                x1      = rnorm(240))
d$outcome <- d$x1 + rnorm(240, sd = ifelse(d$treated == 1, 2, 1))

# Complete randomization: retain the total number treated.
riposte_test(outcome ~ treated, data = d,
             statistic = "quadratic", engine = "permute", seed = 1)
riposte_test(outcome ~ treated, data = d,
             statistic = "quadratic", engine = "asymptotic")

# Blocked randomization: retain each block's number treated.
riposte_test(outcome ~ treated | block, data = d,
             statistic = "cauchy", engine = "permute", seed = 1)
riposte_test(outcome ~ treated | block, data = d,
             statistic = "cauchy", engine = "asymptotic")
```

To test with one representation alone, pass it as the only member of
`representations`. With the raw outcome, every choice of `statistic`
gives the same two-sided permutation test of the treated units’ outcome
sum, which has the same p-value as the difference in means within
blocks:

``` r

riposte_test(outcome ~ treated | block, data = d,
             representations = list(raw = function(y) y), seed = 1)

# one-sided: are treated outcomes higher?
riposte_test(outcome ~ treated | block, data = d,
             representations = list(raw = function(y) y),
             statistic = "max", alternative = "greater", seed = 1)
```

On these data, where treatment changes the spread and not the mean, the
raw outcome alone gives p = 0.53 two-sided, while the blocked Cauchy
combination of all six representations above gives p = 0.004.

Both engines support both designs and both combinations shown above. The
asymptotic engine uses no random draws and ignores `nresample` and
`seed`. It requires `statistic = "quadratic"` or `"cauchy"` and no
covariance adjustment. The default screen, the max, and min-p remain
available with the permutation engine. `coin` is used for comparison
tests, not required to run these approximations.
[`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
continues to report permutation mid-p values, even when used alongside
an asymptotic combined test.

## Status

Version 0.0.0.9004 adds the asymptotic option and unblocked formulas.
`devtools::check()` completed on 2026-09-29 with 0 errors, 0 warnings,
and 0 notes. All 560 test assertions passed, including the size and
power simulations, with no skipped tests (R 4.6.1 on macOS).
Implemented: the outcome representations (open and user-supplied), the
closed-form Strasser-Weber moments (verified against `coin` and against
full enumeration), the quadratic/Cauchy/max combinations, the screen,
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
and
[`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md),
controls-only covariance adjustment (refit per permutation),
cluster-randomized designs, and the combined Stephenson rank test. Size
and power directions are confirmed by simulation, exactness by
enumeration, and correctness by two rounds of adversarial multi-agent
audit.

Open work: the adjusted quadratic/screen use Monte-Carlo moments (the
closed-form moments do not apply under per-draw residualization);
`propertee` interoperation and the few-cluster Neyman SE are not yet
implemented. See `NEWS.md`.

``` r

riposte_test(outcome ~ treated | block, data = d)                  # combined test (screen default)
riposte_test(outcome ~ treated | block, data = d, adjust = ~ x1)   # covariance-adjusted
riposte_test(outcome ~ treated | block, data = d, clusters = "cl") # cluster-randomized
riposte_components(outcome ~ treated | block, data = d)            # per-representation diagnostics
```

See
[`vignette("riposte")`](https://bowers-illinois-edu.github.io/riposte/articles/riposte.md)
for a worked, policy-oriented introduction.

## Related packages

- **manytestsr** — tree-structured testing that locates effects. riposte
  does not do tree testing.
- **CMRSS** — the combined Stephenson rank test.
- **propertee** — design-based estimation. riposte interoperates with it
  for covariance adjustment rather than duplicating it.

## License

MIT (c) Jake Bowers and Myla Burton.
