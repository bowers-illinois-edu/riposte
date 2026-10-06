## test-adjust-validity.R --- the exactness and power claims for adjustment.
##
## THE TWO CLAIMS that justify refit-per-permutation (spec section 8):
##   1. EXACTNESS FOR ANY LEARNER. Refitting the learner inside every
##      re-randomization holds the 0.05 level under the sharp null, even for a
##      learner that overfits (OLS with many covariates). The fit-once /
##      fixed-residual shortcut does NOT: it is anti-conservative when the learner
##      overfits, because the residuals then carry information about the observed
##      assignment. This is the regression test that keeps the shortcut out.
##   2. POWER. When covariates predict the outcome, adjustment increases power to
##      detect a real effect relative to no adjustment.
##
## Monte-Carlo, so slow and skipped on CRAN. Calibrated 2026-06-23 (see the rates
## next to each assertion); thresholds sit several Monte-Carlo SEs from the
## calibrated values.

assign_within_block <- function(block) {
  z <- integer(length(block))
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(length(ix), length(ix) %/% 2L)] <- 1L
  }
  z
}

adj_size <- function(nsims, refit, B = 6L, nb = 10L, d = 20L, nres = 99L) {
  block <- factor(rep(seq_len(B), each = nb)); N <- B * nb
  raw1 <- list(raw = function(y) y)
  fml <- stats::as.formula(paste0("~", paste0("x", seq_len(d), collapse = "+")))
  p <- numeric(nsims)
  for (i in seq_len(nsims)) {
    X <- matrix(rnorm(N * d), N, d)
    y <- rnorm(N)                                    # sharp null: no effect, X unrelated
    z <- assign_within_block(block)
    dat <- data.frame(Y = y, trt = z, blk = block, X)
    names(dat)[-(1:3)] <- paste0("x", seq_len(d))
    p[i] <- riposte_test(Y ~ trt | blk, dat, statistic = "max", representations = raw1,
                         adjust = fml, learner = riposte_lm_learner(),
                         adjust_refit = refit, nresample = nres)$p.value
  }
  mean(p <= 0.05)
}

test_that("refit-per-permutation holds size with an overfitting learner; fit-once does not", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260623)
  ## Calibrated (n = 120, OLS with d = 20 on n_control = 30):
  ##   refit size 0.025, fit-once size 0.183.
  size_refit <- adj_size(100L, refit = TRUE)
  size_fixed <- adj_size(100L, refit = FALSE)
  expect_lt(size_refit, 0.10)                        # refit holds the level (calibrated 0.025)
  expect_gt(size_fixed, 0.12)                        # fit-once inflates it (calibrated 0.183)
  expect_gt(size_fixed - size_refit, 0.05)           # and clearly worse than refit
})

adj_power <- function(nsims, adjusted, B = 8L, nb = 10L, d = 3L, tau = 0.5, nres = 99L) {
  block <- factor(rep(seq_len(B), each = nb)); N <- B * nb
  raw1 <- list(raw = function(y) y)
  fml <- stats::as.formula(paste0("~", paste0("x", seq_len(d), collapse = "+")))
  p <- numeric(nsims)
  for (i in seq_len(nsims)) {
    X <- matrix(rnorm(N * d), N, d)
    z <- assign_within_block(block)
    y <- as.numeric(X %*% c(2, -1.5, 1)) + rnorm(N, sd = 0.7) + tau * z   # covariates predictive
    dat <- data.frame(Y = y, trt = z, blk = block, X)
    names(dat)[-(1:3)] <- paste0("x", seq_len(d))
    p[i] <- if (adjusted)
      riposte_test(Y ~ trt | blk, dat, statistic = "max", representations = raw1,
                   adjust = fml, learner = riposte_ridge_learner(), nresample = nres)$p.value
    else
      riposte_test(Y ~ trt | blk, dat, statistic = "max", representations = raw1,
                   nresample = nres)$p.value
  }
  mean(p <= 0.05)
}

adj_size_combo <- function(nsims, statistic, B = 6L, nb = 10L, d = 6L, nres = 99L) {
  block <- factor(rep(seq_len(B), each = nb)); N <- B * nb
  fml <- stats::as.formula(paste0("~", paste0("x", seq_len(d), collapse = "+")))
  p <- numeric(nsims)
  for (i in seq_len(nsims)) {
    X <- matrix(rnorm(N * d), N, d)
    y <- rnorm(N)                                    # sharp null, X unrelated
    z <- assign_within_block(block)
    dat <- data.frame(Y = y, trt = z, blk = block, X)
    names(dat)[-(1:3)] <- paste0("x", seq_len(d))
    p[i] <- riposte_test(Y ~ trt | blk, dat, statistic = statistic,
                         adjust = fml, learner = riposte_ridge_learner(),
                         nresample = nres)$p.value     # default six representations
  }
  mean(p <= 0.05)
}

test_that("the adjusted quadratic, max, and screen hold the level under the null", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260625)
  ## With POOLED moments (observed + draws) the adjusted combinations are
  ## exactly level-valid -- the metric is a symmetric function of the exchangeable
  ## set, so the combined statistic is permutation-equivariant (see
  ## dev/theory-adjusted-exactness.md). nresample is deliberately small (49), where
  ## the old draws-only metric was badly anti-conservative (~0.12); pooled holds
  ## ~0.05. Two-sided band ~ +/- 3.5 Monte-Carlo SE at n = 200.
  for (st in c("quadratic", "max", "screen")) {
    rate <- adj_size_combo(200L, st, nres = 49L)
    expect_lt(rate, 0.10)
    expect_gt(rate, 0.015)
  }
})

test_that("the adjusted min-p holds the level under the null", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261005)
  ## Min-p needs more re-randomizations than the other combinations to reach
  ## 0.05 at all: each representation's most extreme re-randomization ties for
  ## the smallest mid-p value, so with K representations min-p's p-value cannot
  ## fall below about K / (nresample + 1). At nresample = 49 that floor is 0.10
  ## for the default six, and min-p never rejects; at 199 it is 0.03.
  rate <- adj_size_combo(200L, "minp", nres = 199L)
  expect_lt(rate, 0.10)
  expect_gt(rate, 0.015)
})

test_that("the adjusted quadratic is exact under full enumeration (pooled moments)", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260626)
  ## 2 blocks of 4, 2 treated each -> 36 assignments. Refit per assignment, pooled
  ## quadratic; under H0 the p-value's rank is uniform, so P(p <= k/36) = k/36.
  block <- factor(rep(1:2, each = 4))
  per <- combn(4, 2, function(t) { z <- integer(4); z[t] <- 1L; z }, simplify = FALSE)
  gr <- expand.grid(a = 1:6, b = 1:6)
  G <- vapply(seq_len(nrow(gr)), function(r) c(per[[gr$a[r]]], per[[gr$b[r]]]), numeric(8))
  reps <- riposte_reps_default()
  nsims <- 1500L; p <- numeric(nsims)
  for (i in seq_len(nsims)) {
    X <- matrix(rnorm(8 * 3), 8, 3); y <- rnorm(8)        # sharp null, X unrelated
    obs <- sample.int(ncol(G), 1)
    Tmat <- riposte:::riposte_residual_perm_stats(y, X, block, G,
              riposte_ridge_learner(lambda = 1), reps, refit = TRUE)$stats
    mu <- colMeans(Tmat); S <- stats::cov(Tmat)
    cen <- sweep(Tmat, 2, mu); q <- rowSums((cen %*% MASS::ginv(S)) * cen)
    p[i] <- mean(q >= q[obs])
  }
  ## the empirical CDF should track the uniform-on-grid CDF (SE ~ 0.012 at n=1500)
  expect_lt(abs(mean(p <= 6 / 36 + 1e-9) - 6 / 36), 0.04)
  expect_lt(abs(mean(p <= 18 / 36 + 1e-9) - 18 / 36), 0.05)
  expect_lt(abs(mean(p) - (ncol(G) + 1) / (2 * ncol(G))), 0.04)
})

test_that("adjustment with predictive covariates increases power", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260624)
  ## Calibrated (n = 120, tau = 0.5): unadjusted 0.12, adjusted 0.82.
  pow_unadj <- adj_power(100L, adjusted = FALSE)
  pow_adj <- adj_power(100L, adjusted = TRUE)
  expect_gt(pow_adj, 0.50)
  expect_gt(pow_adj - pow_unadj, 0.20)
})
