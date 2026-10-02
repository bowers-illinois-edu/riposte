## test-pvalues.R
##
## The statistical content: the mid-p value must be strictly inside (0, 1) so the
## Cauchy transform never reaches a pole (the bug the spec calls out --- a plain
## permutation p-value can hit exactly 1, which a naive clamp turns into a huge
## negative term that swamps the average). The ACAT term must reproduce the plain
## tan transform away from the poles and stay finite at them.

test_that("mid-p values are strictly inside (0, 1)", {
  set.seed(1)
  v <- rnorm(200)
  p <- riposte_midp(v)
  expect_true(all(p > 0 & p < 1))
  ## extremes: the most extreme |v| gets 0.5/n, the least gets (n - 0.5)/n
  n <- length(v)
  expect_equal(min(p), 0.5 / n)
  expect_equal(max(p), (n - 0.5) / n)
})

test_that("mid-p ranks by |v|: larger magnitude -> smaller p", {
  v <- c(0.1, -3, 2, -0.5)
  p <- riposte_midp(v)
  expect_equal(order(p), order(abs(v), decreasing = TRUE))
})

test_that("mid-p keeps the Cauchy transform finite where a plain p-value would not", {
  ## a statistic vector whose least extreme entry would get a plain two-sided
  ## p-value of exactly 1; mid-p keeps it inside (0, 1) so tan() is finite
  v <- c(5, 4, 3, 0)
  p <- riposte_midp(v)
  expect_true(all(is.finite(tan((0.5 - p) * pi))))
})

test_that("the Cauchy combination survives a representation sitting at the plain-p pole", {
  ## column 1's observed value (row 1) is the LEAST extreme (|T| = 0), so its plain
  ## two-sided permutation p-value is exactly 1 -- the pole where tan() diverges.
  ## The mid-p the Cauchy uses keeps it strictly inside (0, 1), so the combination
  ## stays finite where a clamp-on-the-plain-p combiner would blow up.
  Tmat <- cbind(c(0, 1, -2, 3, -1, 2), c(1.5, -0.3, 0.7, -1.1, 0.2, -0.8))
  res <- riposte:::riposte_cauchy_from_T(Tmat)
  expect_true(is.finite(res$statistic))
  expect_true(res$p.value > 0 && res$p.value <= 1)

  plain_p <- mean(abs(Tmat[, 1]) >= abs(Tmat[1, 1]))   # the observed is least extreme
  expect_equal(plain_p, 1)                             # exactly at the pole
  clamp_term <- tan((0.5 - pmin(plain_p, 1 - 1e-15)) * pi)
  expect_gt(abs(clamp_term), 1e10)                     # a clamp would diverge here

  midp_obs <- riposte_midp(Tmat[, 1])[1]               # mid-p stays interior
  expect_true(midp_obs > 0 && midp_obs < 1)
  expect_lt(abs(riposte_acat_term(midp_obs)), 10)      # and its transform is bounded
})

test_that("acat_term equals tan((0.5 - p) pi) away from the poles", {
  p <- c(0.2, 0.5, 0.8)
  expect_equal(riposte_acat_term(p), tan((0.5 - p) * pi))
})

test_that("acat_term stays finite at the poles via the small-/large-p limits", {
  expect_true(is.finite(riposte_acat_term(1e-300)))
  expect_true(is.finite(riposte_acat_term(1 - 1e-16)))
  ## sign is right: tiny p -> large positive, p near 1 -> large negative
  expect_gt(riposte_acat_term(1e-8), 0)
  expect_lt(riposte_acat_term(1 - 1e-8), 0)
})

test_that("values equal up to floating-point error count as ties", {
  ## Rank scores give many re-randomizations whose treated-score sums are
  ## exactly equal in exact arithmetic, but floating-point sums of the same
  ## numbers in a different order can differ in the last bits. 0.1 + 0.2 and 0.3
  ## are the textbook pair. A tied re-randomization must count as at least as
  ## extreme as the observed one; if last-bit noise put it just below, the
  ## p-value would be too small and the test anti-conservative.
  tied <- 0.1 + 0.2                               # 0.30000000000000004
  stat <- c(tied, 0.3, 0.3, 0.1)
  expect_equal(riposte:::riposte_perm_pvalue(stat), 3 / 4)

  ## mid-p: the three tied values share one value, each with half the tied mass
  p <- riposte_midp(c(0.3, tied, 0.3, 1), "greater")
  expect_equal(p[1], p[2])
  expect_equal(p[1], (1 + 0.5 * 3) / 4)
  ## values farther apart than floating-point error are not ties
  expect_lt(riposte_midp(c(0.3, 0.3 + 1e-6), "greater")[2],
            riposte_midp(c(0.3, 0.3 + 1e-6), "greater")[1])
})
