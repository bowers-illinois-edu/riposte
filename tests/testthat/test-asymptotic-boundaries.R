## Numerical boundaries matter for the analytic Cauchy calculation: an
## individual p-value rounded to 1 gives an infinite negative tangent, while
## a tiny p-value can give a tangent too large to store as a double.
## Approximation mode must not silently switch to random permutations.

test_that("an exactly zero statistic has asymptotic p-value one", {
  d <- data.frame(Y = c(1, 2, 3, 4), trt = c(1, 0, 0, 1))
  set.seed(12)
  before <- .Random.seed
  for (st in c("quadratic", "cauchy")) {
    res <- riposte_test(Y ~ trt, d, statistic = st,
                        representations = list(raw = identity),
                        engine = "asymptotic")
    expect_equal(res$p.value, 1)
    expect_identical(.Random.seed, before)
  }
})

test_that("a p-value near one retains its finite Cauchy statistic", {
  ## The -1 and 1 sit in the control group so the treated sum has no
  ## cancellation. With both in the treated group, the sum depended on BLAS
  ## summation order: OpenBLAS (Ubuntu CI) adds (-1 + -1e-20) + 1 = 0, where
  ## reference BLAS adds (-1 + 1) + -1e-20 = -1e-20. A zero sum is the exact
  ## p = 1 case, whose statistic is -Inf by design, not the case tested here.
  d <- data.frame(Y = c(-1e-20, 1e-20, 1, -1), trt = c(1, 0, 0, 0))
  ## The mean is zero, the treated sum is -1e-20, and its permutation
  ## variance is 1 * 3 / (4 * 3) * 2 = 0.5, so the squared standardized sum
  ## is x = (1e-20)^2 / 0.5 = 2e-40. Compute the small lower chi-square tail
  ## directly: subtracting the upper tail from 1 loses it. The expected value
  ## does not call pchisq(). Where pchisq() is inaccurate at so small an x,
  ## it would make the expected value wrong in the same way as the code.
  ## P(chi2_1 < x) = P(|Z| < sqrt(x)) = sqrt(2x/pi) * (1 - x/6 + ...), and at
  ## x = 2e-40 the factor (1 - x/6) equals 1 in double precision.
  x <- 1e-40 / 0.5
  lower <- sqrt(2 * x / pi)
  expected <- -1 / tan(pi * lower)
  ## this checks the precision of Liu and Xie's untruncated conversion near
  ## p = 1, so it asks for that conversion with cauchy_truncation = 1
  res <- riposte_test(Y ~ trt, d, statistic = "cauchy",
                      representations = list(raw = identity),
                      engine = "asymptotic", cauchy_truncation = 1)
  expect_true(is.finite(res$statistic))
  expect_equal(res$statistic / expected, 1, tolerance = 1e-12)
})

## Below x = 1e-30, riposte_log_chisq1_lower() computes log P(chi2_1 < x)
## from a formula instead of calling pchisq(). Myla Burton reported
## platform-dependent results from pchisq() at extreme small x; we have not
## reproduced them. The formula comes from the series above with the factor
## (1 - x/6) dropped, which equals 1 in double precision once x is below
## about 1e-16.

test_that("the small-x formula agrees with pchisq where pchisq is accurate", {
  ## tol = 1 makes the function use the formula at values where it would
  ## otherwise call pchisq(), so this checks the formula itself against an
  ## independent calculation.
  x <- c(1e-28, 1e-24, 1e-20)
  expect_equal(riposte_log_chisq1_lower(x, tol = 1),
               pchisq(x, df = 1, log.p = TRUE), tolerance = 1e-12)
})

test_that("Cauchy preserves a tiny tail even if its statistic overflows", {
  d <- data.frame(Y = rep(0:1, each = 712), trt = rep(0:1, each = 712))
  expected <- pchisq(1423, df = 1, lower.tail = FALSE)
  expect_gt(expected, 0)
  expect_lt(expected, .Machine$double.xmin)
  ## Identical component p-values must combine to that same p-value.
  res <- riposte_test(Y ~ trt, d, statistic = "cauchy",
                      representations = list(raw = identity,
                                             twice = function(y) 2 * y),
                      engine = "asymptotic")
  expect_gt(res$p.value, 0)
  expect_equal(res$p.value / expected, 1, tolerance = 1e-9)
})

test_that("scores varying only in nonrandomized blocks add no information", {
  d <- data.frame(Y = c(1, 2, 4, 8, 11, 12, 14, 18),
                  trt = c(1, 0, 0, 1, 1, 1, 1, 1),
                  blk = factor(rep(1:2, each = 4)))
  reps <- list(raw = identity,
               fixed = function(y) if (min(y) < 10) rep(0, length(y)) else y)
  for (st in c("quadratic", "cauchy")) {
    raw <- riposte_test(Y ~ trt | blk, d, statistic = st,
                        representations = list(raw = identity),
                        engine = "asymptotic")
    res <- riposte_test(Y ~ trt | blk, d, statistic = st,
                        representations = reps, engine = "asymptotic")
    expect_equal(res$kept, "raw")
    expect_equal(res$dropped, "fixed")
    expect_equal(res$p.value, raw$p.value)
  }
})

test_that("asymptotic mode ignores the unused resampling count and seed", {
  d <- asymptotic_example()
  set.seed(19)
  before <- .Random.seed
  res <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                      engine = "asymptotic", nresample = 0, seed = 99)
  expect_equal(res$nresample, 0L)
  expect_identical(.Random.seed, before)
})

test_that("unsupported approximations give an explicit error", {
  d <- asymptotic_example()
  ## A maximum needs a joint reference distribution. A shrunk quadratic
  ## does not generally have the ordinary chi-square reference. Neither is
  ## part of the requested unadjusted quadratic and Cauchy implementation.
  for (st in c("max", "screen")) {
    expect_error(riposte_test(Y ~ trt | blk, d, statistic = st,
                              engine = "asymptotic"),
                 'asymptotic.*quadratic.*cauchy')
  }
  ## A model refit under each assignment needs its own approximation theory.
  expect_error(riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                            engine = "asymptotic", adjust = ~ Y),
               "asymptotic.*adjustment")
})
