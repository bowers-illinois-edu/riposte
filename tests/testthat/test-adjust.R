## test-adjust.R --- covariance adjustment plumbing (fast unit tests).
##
## The statistical content here is the contract, not the size/power claims (those
## are the slow simulations in test-adjust-validity.R):
##   - the ridge learner is controls-only (never sees treated outcomes) and equals
##     the closed-form ridge solution;
##   - refit-per-permutation and the fit-once shortcut agree on the OBSERVED
##     assignment (both fit on the observed controls) but differ on the draws;
##   - riposte_test wires the adjustment through and reports it.

make_xy <- function(seed = 1, N = 40, d = 3) {
  set.seed(seed)
  X <- matrix(rnorm(N * d), N, d)
  beta <- rnorm(d)
  list(X = X, Y = as.numeric(X %*% beta) + rnorm(N), control = rep(c(TRUE, FALSE), length.out = N))
}

test_that("ridge learner equals the closed-form controls-only ridge", {
  dat <- make_xy()
  X <- dat$X; Y <- dat$Y; ctrl <- dat$control; lam <- 2.5
  yhat <- riposte_ridge_learner(lambda = lam)(X, Y, ctrl)
  ## closed form: centre on control means, ridge without intercept, add it back
  Xc <- X[ctrl, , drop = FALSE]; Yc <- Y[ctrl]
  xm <- colMeans(Xc); ym <- mean(Yc)
  beta <- solve(crossprod(sweep(Xc, 2, xm)) + lam * diag(ncol(X)),
                crossprod(sweep(Xc, 2, xm), Yc - ym))
  ref <- as.numeric(sweep(X, 2, xm) %*% beta) + ym
  expect_equal(yhat, ref)
})

test_that("the learner is controls-only: treated outcomes do not affect predictions", {
  dat <- make_xy()
  yhat1 <- riposte_ridge_learner(lambda = 1)(dat$X, dat$Y, dat$control)
  Y2 <- dat$Y
  Y2[!dat$control] <- Y2[!dat$control] + 100        # corrupt treated outcomes
  yhat2 <- riposte_ridge_learner(lambda = 1)(dat$X, Y2, dat$control)
  expect_equal(yhat1, yhat2)
})

test_that("ridge CV picks a finite penalty and predicts sensibly", {
  dat <- make_xy(N = 60, d = 4)
  yhat <- riposte_ridge_learner(lambda = NULL)(dat$X, dat$Y, dat$control)
  expect_length(yhat, nrow(dat$X))
  expect_true(all(is.finite(yhat)))
})

test_that("ridge LOOCV includes the intercept leverage and does not collapse to the smallest penalty", {
  ## Regression test for the audit's major finding: the SVD LOOCV must carry the
  ## unpenalized intercept's 1/n_control leverage. Without it, CV picks the
  ## smallest penalty in the high-dimensional regime (d ~ n_control), where ridge
  ## should regularize hardest. Check against honest brute-force LOOCV that refits
  ## and re-centres on each leave-one-out training set.
  set.seed(2024)
  nc <- 20L; d <- 19L                          # d = n_control - 1: near-interpolation
  Xc <- matrix(rnorm(nc * d), nc, d); Yc <- rnorm(nc)   # pure noise: nothing to fit
  grid <- 10^seq(-3, 3, length.out = 50)
  honest <- vapply(grid, function(lam) {
    err <- numeric(nc)
    for (i in seq_len(nc)) {
      Xtr <- Xc[-i, , drop = FALSE]; Ytr <- Yc[-i]
      xm <- colMeans(Xtr); ym <- mean(Ytr)
      Xcc <- sweep(Xtr, 2, xm); Ycc <- Ytr - ym
      b <- solve(crossprod(Xcc) + lam * diag(d), crossprod(Xcc, Ycc))
      err[i] <- (Yc[i] - (as.numeric((Xc[i, ] - xm) %*% b) + ym))^2
    }
    mean(err)
  }, numeric(1))
  lam_honest <- grid[which.min(honest)]
  ## the package CV is called on globally-centred control data (as the learner does)
  lam_pkg <- riposte:::riposte_ridge_cv(sweep(Xc, 2, colMeans(Xc)), Yc - mean(Yc), grid)
  expect_equal(lam_pkg, lam_honest)
  expect_gt(lam_pkg, min(grid))                # not collapsed to the smallest penalty
})

test_that("refit and fit-once agree on the observed row but differ on the draws", {
  set.seed(3)
  block <- factor(rep(1:4, each = 8)); N <- length(block)
  z <- integer(N); for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(8, 4)] <- 1L }
  X <- matrix(rnorm(N * 3), N, 3)
  y <- as.numeric(X %*% c(1, -1, 0.5)) + rnorm(N)
  G <- riposte_block_draws(z, block, nresample = 30)
  reps <- list(raw = function(v) v)
  refit <- riposte:::riposte_residual_perm_stats(y, X, block, G, riposte_ridge_learner(1), reps, refit = TRUE)
  fixed <- riposte:::riposte_residual_perm_stats(y, X, block, G, riposte_ridge_learner(1), reps, refit = FALSE)
  ## observed assignment (row 1): both fit on the observed controls -> identical
  expect_equal(refit$stats[1, ], fixed$stats[1, ])
  ## the draws use different control sets under refit, so the statistics differ
  expect_false(isTRUE(all.equal(refit$stats[-1, ], fixed$stats[-1, ])))
})

test_that("the adjusted refit engine selects representations symmetrically", {
  ## Regression test for the exactness gap the adversarial review found: the kept
  ## representation set must be a function of the whole orbit, not of the observed
  ## assignment alone, or the combined statistic depends on which assignment is
  ## "observed" and exactness breaks. Here a small design where a representation can
  ## be within-block-constant on some assignments' residuals; the kept set must not
  ## change when a different assignment plays the role of observed (column 1).
  set.seed(3)
  block <- factor(rep(1:2, each = 4))
  per <- combn(4, 2, function(t) { z <- integer(4); z[t] <- 1L; z }, simplify = FALSE)
  gr <- expand.grid(a = 1:6, b = 1:6)
  G <- vapply(seq_len(nrow(gr)), function(r) c(per[[gr$a[r]]], per[[gr$b[r]]]), numeric(8))
  X <- matrix(rnorm(8 * 3), 8, 3); y <- rnorm(8); reps <- riposte_reps_default()
  eng <- riposte:::riposte_residual_perm_stats
  k1 <- eng(y, X, block, G[, c(1, 2:36)], riposte_ridge_learner(1), reps)$kept
  k2 <- eng(y, X, block, G[, c(5, (1:36)[-5])], riposte_ridge_learner(1), reps)$kept
  expect_setequal(k1, k2)
  ## and the pooled covariance of the kept columns has no zero-variance column,
  ## so the screen's cov2cor cannot hit a divide-by-zero (the second finding)
  Tmat <- eng(y, X, block, G, riposte_ridge_learner(1), reps)$stats
  expect_true(all(apply(Tmat, 2, stats::var) > 1e-12))
})

test_that("riposte_test wires adjustment through and reports it", {
  set.seed(4)
  block <- factor(rep(1:5, each = 8)); N <- length(block)
  z <- integer(N); for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(8, 4)] <- 1L }
  X <- matrix(rnorm(N * 2), N, 2)
  y <- as.numeric(X %*% c(1.5, -1)) + rnorm(N) + 0.4 * z
  d <- data.frame(Y = y, trt = z, blk = block, x1 = X[, 1], x2 = X[, 2])
  r <- riposte_test(Y ~ trt | blk, d, statistic = "max",
                    adjust = ~ x1 + x2, nresample = 199, seed = 9)
  expect_s3_class(r, "riposte_test")
  expect_true(isTRUE(r$adjusted))
  expect_true(r$p.value > 0 && r$p.value <= 1)
  expect_output(print(r), "covariance-adjusted")
})

test_that("a user-supplied learner flows through", {
  set.seed(5)
  block <- factor(rep(1:4, each = 8)); N <- length(block)
  z <- integer(N); for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(8, 4)] <- 1L }
  d <- data.frame(Y = rnorm(N), trt = z, blk = block, x1 = rnorm(N))
  ## a trivial learner: predict the control mean for everyone
  mean_learner <- function(X, Y, control) rep(mean(Y[control]), length(Y))
  r <- riposte_test(Y ~ trt | blk, d, statistic = "max",
                    adjust = ~ x1, learner = mean_learner, nresample = 99, seed = 1)
  expect_s3_class(r, "riposte_test")
  expect_true(r$p.value > 0 && r$p.value <= 1)
})

test_that("the covariate matrix drops the intercept", {
  d <- data.frame(x1 = rnorm(5), x2 = rnorm(5))
  M <- riposte:::riposte_covariate_matrix(~ x1 + x2, d)
  expect_false("(Intercept)" %in% colnames(M))
  expect_setequal(colnames(M), c("x1", "x2"))
})
