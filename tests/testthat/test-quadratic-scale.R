## test-quadratic-scale.R
##
## Issue #1: the quadratic statistic Q = (T - mu)' Sigma^+ (T - mu) must not
## depend on the units of the outcome, and its df must count the directions Q
## actually uses.
##
## The statistical content:
##   - A test of the sharp null is a statement about the outcomes and the
##     assignment, not about the units the outcomes are recorded in. Multiplying
##     Y by a constant leaves the rank scores unchanged and rescales the raw and
##     distance scores, so Q, its df, and its p-value must not move.
##   - For an outcome with K <= 7 distinct values in one block, the six default
##     representations are functions of K values, so their centered score
##     vectors span exactly K - 1 directions and fill them. Q is then the
##     quadratic form of the treated counts at each value, which is Pearson's
##     chi-square statistic for the 2 x K table times (n - 1)/n, with K - 1 df.
##     Inverting the raw covariance with MASS::ginv() fails this on survey
##     scales: the rank sum's variance is ~10^7 times the raw sum's, and ginv()
##     drops a real direction whose eigenvalue falls below its relative cutoff.
##   - Where ginv() on the raw covariance already keeps every direction, the
##     fix must give the same Q as before (the change is numerical only).

## a 5-point outcome on 0, 0.25, ..., 1 with 962 units: the case in issue #1,
## where ginv() on the raw covariance keeps 3 of the 4 directions
five_point <- function() {
  set.seed(1)
  n <- 962
  y <- sample(c(0, 0.25, 0.5, 0.75, 1), n, replace = TRUE, prob = c(.34, .14, .19, .17, .16))
  z <- sample(rep(0:1, c(310, 652)))
  data.frame(Y = y, Z = z)
}

pearson_adjusted <- function(d) {
  n <- nrow(d)
  unname(suppressWarnings(stats::chisq.test(table(d$Z, d$Y), correct = FALSE))$statistic) * (n - 1) / n
}

test_that("the permutation quadratic does not depend on the units of the outcome", {
  d <- five_point()
  fits <- lapply(c(1, 100, 0.01), function(k)
    riposte_test(Y ~ Z, data = transform(d, Y = Y * k), statistic = "quadratic", seed = 1))
  for (f in fits[-1]) {
    expect_equal(f$statistic, fits[[1]]$statistic, tolerance = 1e-8)
    expect_equal(f$df, fits[[1]]$df)
    expect_equal(f$p.value, fits[[1]]$p.value)
  }
})

test_that("on a K-valued outcome the quadratic is Pearson's chi-square times (n - 1)/n on K - 1 df", {
  d <- five_point()
  f <- riposte_test(Y ~ Z, data = d, statistic = "quadratic", seed = 1)
  expect_equal(f$statistic, pearson_adjusted(d), tolerance = 1e-8)
  expect_equal(f$df, 4L)
  ## a 4-point outcome: 3 directions
  set.seed(2)
  d4 <- data.frame(Y = sample(0:3 / 3, 600, replace = TRUE, prob = c(.25, .2, .37, .18)),
                   Z = sample(rep(0:1, 300)))
  f4 <- riposte_test(Y ~ Z, data = d4, statistic = "quadratic", seed = 1)
  expect_equal(f4$statistic, pearson_adjusted(d4), tolerance = 1e-8)
  expect_equal(f4$df, 3L)
})

test_that("the asymptotic quadratic does not depend on the units and matches Pearson with K - 1 df", {
  d <- five_point()
  fits <- lapply(c(1, 100, 0.01), function(k)
    riposte_test(Y ~ Z, data = transform(d, Y = Y * k), statistic = "quadratic", engine = "asymptotic"))
  for (f in fits) {
    expect_equal(f$statistic, pearson_adjusted(d), tolerance = 1e-8)
    expect_equal(f$df, 4L)
    expect_equal(f$p.value, stats::pchisq(pearson_adjusted(d), 4, lower.tail = FALSE), tolerance = 1e-8)
  }
})

test_that("rescaling the outcome leaves a blocked quadratic unchanged", {
  set.seed(3)
  block <- factor(rep(1:4, each = 150))
  y <- sample(1:5, 600, replace = TRUE)
  z <- unlist(lapply(split(seq_along(y), block), function(ix) sample(rep(0:1, length(ix) / 2))))
  d <- data.frame(Y = y, Z = z, B = block)
  a <- riposte_test(Y ~ Z | B, data = d, statistic = "quadratic", seed = 1)
  b <- riposte_test(Y ~ Z | B, data = transform(d, Y = Y / 1000), statistic = "quadratic", seed = 1)
  expect_equal(b$statistic, a$statistic, tolerance = 1e-8)
  expect_equal(b$df, a$df)
  expect_equal(b$p.value, a$p.value)
})

test_that("where ginv() on the raw covariance kept every direction, Q is unchanged", {
  ## a continuous outcome in six blocks: no eigenvalue of Sigma is anywhere
  ## near ginv()'s cutoff, so the old and new forms are the same number
  set.seed(20260623)
  block <- factor(rep(1:6, each = 10))
  z <- unlist(lapply(split(seq_along(block), block), function(ix) sample(rep(0:1, 5))))
  y <- rnorm(60) + 0.5 * z
  sm <- riposte_score_matrix(y, block)
  mom <- riposte_sw_moments(sm$scores, z, block)
  cen <- as.numeric(crossprod(z, sm$scores)) - mom$mu
  q_old <- drop(t(cen) %*% MASS::ginv(mom$Sigma) %*% cen)
  f <- riposte_test(Y ~ Z | B, data = data.frame(Y = y, Z = z, B = block), statistic = "quadratic", seed = 1)
  expect_equal(f$statistic, q_old, tolerance = 1e-8)
  expect_equal(f$df, qr(mom$Sigma)$rank)
})
