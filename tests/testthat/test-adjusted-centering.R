## test-adjusted-centering.R
##
## THE POINT. With covariance adjustment riposte refits the residuals at every
## re-randomization, so the null mean of each treated sum is its average over
## the observed and redrawn assignments, which is not zero. Jake decided on
## 2026-11-16 that every combination measures departures from that average: the
## max, the quadratic, and the min-p already did, and the Cauchy combination and
## the hybrid now do too. The test of that is translation invariance: adding a
## constant to every treated sum of a representation, at every assignment, moves
## the average by the same constant and must leave every p-value unchanged.

test_that("every adjusted combination is unchanged when a constant is added to every treated sum", {
  set.seed(20261116)
  block <- factor(rep(1:4, each = 8))
  z <- integer(32)
  for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(8, 4)] <- 1L }
  y <- rnorm(32) + 0.6 * z
  X <- cbind(x = y + rnorm(32))
  RNGkind("L'Ecuyer-CMRG"); set.seed(1)
  G <- riposte_block_draws(z, block, 199L)
  rp <- riposte:::riposte_residual_perm_stats(y, X, block, G, riposte_ridge_learner(),
                                              riposte_reps_default(), TRUE)
  Tmat <- rp$stats
  shifted <- sweep(Tmat, 2, seq_len(ncol(Tmat)) * 3.1, "+")
  for (stat in c("quadratic", "max", "minp", "cauchy", "hybrid", "screen")) {
    a <- riposte:::riposte_adjusted_combination(Tmat, stat, riposte_screen_control())
    b <- riposte:::riposte_adjusted_combination(shifted, stat, riposte_screen_control())
    expect_equal(b$result$p.value, a$result$p.value, label = stat)
    expect_equal(b$result$statistic, a$result$statistic, tolerance = 1e-8, label = stat)
  }
  ## one-sided Cauchy as well
  a <- riposte:::riposte_adjusted_combination(Tmat, "cauchy", riposte_screen_control(), "less")
  b <- riposte:::riposte_adjusted_combination(shifted, "cauchy", riposte_screen_control(), "less")
  expect_equal(b$result$p.value, a$result$p.value)
})

test_that("without adjustment the Cauchy and hybrid p-values are what they were", {
  ## the null mean of a treated sum of within-block-centered scores is zero, so
  ## measuring from it changes nothing; compare with mid-p values from zero
  set.seed(3)
  block <- factor(rep(1:3, each = 10))
  z <- integer(30)
  for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(10, 5)] <- 1L }
  y <- rnorm(30)
  sm <- riposte_score_matrix(y, block)
  RNGkind("L'Ecuyer-CMRG"); set.seed(2)
  G <- riposte_block_draws(z, block, 199L)
  Tmat <- crossprod(G, sm$scores)
  P0 <- apply(Tmat, 2, riposte_midp)
  expected <- riposte_perm_pvalue(rowMeans(riposte_acat_term(P0)))
  expect_equal(riposte:::riposte_cauchy(sm$scores, z, block, draws = G)$p.value, expected)
})
