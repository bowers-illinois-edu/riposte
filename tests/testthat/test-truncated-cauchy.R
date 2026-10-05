## test-truncated-cauchy.R
##
## THE POINT. Liu and Xie's Cauchy combination converts each p-value p to
## tan((0.5 - p) * pi). A p-value of exactly 1 becomes minus infinity, and one
## such value makes the combined p-value 1 whatever the other p-values are. On
## coarse outcomes with few units, a representation's treated sum often equals
## its mean exactly, so its large-sample p-value is exactly 1 (Wallsten and
## Nteta 2016, 36 respondents: the difference in means gives 0.028 and the
## large-sample Cauchy combination gives 1). A p-value just below 1 does almost
## as much damage: 0.999 converts to about -318.
##
## Gui, Jiang, and Wang (2025, Biometrika 112(4), asaf038) convert p to
## tan((0.5 - t * p) * pi) with t = 0.9, the quantile of a Cauchy distribution
## cut off below tan(-0.4 * pi) = -3.08. No converted value can be less than
## -3.08. The cut-off distribution keeps a share t of the Cauchy's probability,
## so its upper tail is the Cauchy's divided by t, and with n p-values whose
## converted values sum to S the combined p-value is min(1, n * P(Cauchy > S) / t).
## riposte uses this for its large-sample Cauchy combination and hybrid, with
## t = 0.9 by default; cauchy_truncation = 1 gives Liu and Xie's original.

test_that("the conversion sends p = 1 to tan(-0.4 pi), not to minus infinity", {
  ## one p-value of 1 among four small ones: the original combination gives 1,
  ## the truncated one does not
  p <- c(0.01, 0.01, 0.01, 0.01, 1)
  expect_equal(riposte_truncated_cauchy(p, truncation = 1), 1)
  S <- 4 * tan((0.5 - 0.9 * 0.01) * pi) + tan(-0.4 * pi)
  expect_equal(riposte_truncated_cauchy(p), 5 * pcauchy(S, lower.tail = FALSE) / 0.9)
  expect_lt(riposte_truncated_cauchy(p), 0.05)
})

test_that("one p-value combines to itself", {
  ## with n = 1, P(Cauchy > tan((0.5 - t p) pi)) = t p, so the combined p-value
  ## is t p / t = p, for any truncation
  for (p in c(1e-200, 0.003, 0.2, 0.7, 1)) {
    expect_equal(riposte_truncated_cauchy(p), p, tolerance = 1e-12)
    expect_equal(riposte_truncated_cauchy(p, truncation = 0.5), p, tolerance = 1e-12)
  }
})

test_that("a tiny p-value neither overflows nor is lost", {
  ## the converted value of 1e-300 is about 1 / (0.9e-300 pi), too large to
  ## store as tan(), so the combination works with its logarithm. One tiny
  ## p-value among n then gives about n times that p-value, Bonferroni's answer.
  res <- riposte_truncated_cauchy(c(1e-300, 0.5, 0.5))
  expect_gt(res, 0)
  expect_equal(res / 3e-300, 1, tolerance = 1e-6)
})

test_that("the truncation must be in (0, 1]", {
  expect_error(riposte_truncated_cauchy(0.5, truncation = 0), "truncation")
  expect_error(riposte_truncated_cauchy(0.5, truncation = 1.1), "truncation")
})

test_that("riposte_test's large-sample Cauchy combination is the truncated combination of its component p-values", {
  d <- wallsten_36()
  res <- riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic")
  expect_equal(res$p.value, riposte_truncated_cauchy(res$component_p))
  ## the max distance's component p-value is exactly 1 here (its treated mean
  ## 31/36 equals its overall mean), and the raw, rank, and Huber p-values are
  ## near 0.03; the original combination gives 1, the truncated one about 0.057
  expect_equal(unname(res$component_p["max_dist"]), 1)
  expect_equal(res$p.value, 0.0574, tolerance = 1e-3)
  old <- riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic", cauchy_truncation = 1)
  expect_equal(old$p.value, 1)
})

test_that("cauchy_truncation changes only the large-sample engine", {
  ## the re-randomization Cauchy already uses mid-p values, which stay below 1,
  ## and its p-value comes from the draws; the argument leaves it unchanged
  d <- wallsten_36()
  a <- riposte_test(Y ~ Z, d, statistic = "cauchy", seed = 1)
  b <- riposte_test(Y ~ Z, d, statistic = "cauchy", seed = 1, cauchy_truncation = 1)
  expect_equal(a$p.value, b$p.value)
})

test_that("under the sharp null the large-sample truncated Cauchy and hybrid keep the 0.05 level on a coarse small outcome", {
  skip_on_cran()
  ## hold the 36 Wallsten and Nteta outcomes fixed, so the treatment changes no
  ## one's outcome, and reassign treatment 1,000 times; SE of the rejection rate
  ## near 0.05 is sqrt(0.05 * 0.95 / 1000) = 0.0069
  d <- wallsten_36()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261004)
  rej <- t(vapply(seq_len(1000), function(i) {
    d$Z <- sample(d$Z)
    c(cauchy = riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic")$p.value,
      hybrid = riposte_test(Y ~ Z, d, statistic = "hybrid", engine = "asymptotic")$p.value) <= 0.05
  }, logical(2)))
  expect_true(all(colMeans(rej) < 0.05 + 3 * sqrt(0.05 * 0.95 / 1000)))
})
