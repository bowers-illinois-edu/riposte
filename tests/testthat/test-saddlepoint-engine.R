## The optional saddlepoint engine computes a combination without re-randomizing.
## Two combinations are draws-free: the UNADJUSTED Cauchy (each representation's
## two-sided permutation p-value from fastperm's scalar saddlepoint, combined by the
## Liu-Xie analytic tail) and the UNADJUSTED quadratic (the Hansen-Bowers omnibus
## Q = (T - mu)' Sigma^{-1} (T - mu), whose within-block permutation tail comes from
## fastperm's multivariate saddlepoint, M2). These tests fix the contract: the
## saddlepoint marginals must match the brute-force permutation marginals, the
## combined p-values must be valid probabilities tracking their permutation
## counterparts, and the engine must refuse the cases it does not yet support. All
## require fastperm to be installed (and M2 needs a fastperm with
## fastperm_spa_quadratic, i.e. the route-b feature).

## a block-randomized data frame with a sharp signal
mk_df <- function(nb = 20L, bs = 6L, m = 3L, tau = 0.4, seed = 1L) {
  set.seed(seed)
  n <- nb * bs
  block <- rep(seq_len(nb), each = bs)
  z <- numeric(n)
  for (ix in split(seq_len(n), block)) z[ix[seq_len(m)]] <- 1
  data.frame(y = rnorm(n) + tau * z, z = z, block = factor(block))
}

test_that("saddlepoint marginals match the brute-force two-sided mid-p", {
  skip_if_not_installed("fastperm")
  ## a representative design (saddlepoint accuracy is O(n^-3/2), so a too-small
  ## design overstates the gap); 30 blocks of 10 is closer to riposte's use
  d <- mk_df(nb = 30L, bs = 10L, m = 5L)
  sm <- riposte_score_matrix(d$y, d$block)
  spa <- riposte_cauchy_spa(sm$scores, d$z, d$block)

  ## brute-force two-sided mid-p of each representation at the observed assignment
  draws <- riposte_block_draws(d$z, d$block, 15000L)
  Tm <- crossprod(draws, sm$scores)
  brute_midp <- apply(Tm, 2, function(col) riposte_midp(col)[1])

  expect_lt(max(abs(spa$component_p - brute_midp)), 0.02)
})

test_that("the saddlepoint Cauchy p-value is valid and tracks the permutation one", {
  skip_if_not_installed("fastperm")
  d <- mk_df()
  p_perm <- riposte_test(y ~ z | block, d, statistic = "cauchy",
                         engine = "permute", nresample = 4999L)$p.value
  p_spa  <- riposte_test(y ~ z | block, d, statistic = "cauchy",
                         engine = "saddlepoint")$p.value

  expect_true(p_spa > 0 && p_spa <= 1)
  ## same calibration target (both reference the within-block permutation null),
  ## so they agree up to the analytic-vs-permutation difference and MC noise
  expect_lt(abs(p_spa - p_perm), 0.08)
})

test_that("the saddlepoint engine computes a draws-free quadratic that tracks permute", {
  skip_if_not_installed("fastperm")
  skip_if_not(exists("fastperm_spa_quadratic", where = asNamespace("fastperm")),
              "fastperm lacks fastperm_spa_quadratic (needs the route-b feature)")
  ## few representations so the tensor Gauss-Hermite grid is feasible: M2 is the
  ## small-r engine (many representations are the unbuilt sparse-grid/QMC regime)
  reps <- riposte_reps_default()[c("raw", "rank")]
  d <- mk_df(nb = 30L, bs = 10L, m = 5L)
  p_perm <- riposte_test(y ~ z | block, d, statistic = "quadratic",
                         representations = reps, engine = "permute",
                         nresample = 4999L)$p.value
  res_spa <- riposte_test(y ~ z | block, d, statistic = "quadratic",
                          representations = reps, engine = "saddlepoint")

  expect_true(res_spa$p.value > 0 && res_spa$p.value <= 1)
  expect_identical(res_spa$combination, "quadratic")
  ## both reference the same within-block permutation null and the SAME statistic
  ## Q; M2 carries the true non-normality of Q and is conservative, so it agrees
  ## with the permutation quadratic up to the few-percent M2 gap and the MC noise
  expect_lt(abs(res_spa$p.value - p_perm), 0.08)
})

test_that("the saddlepoint screen chooses a draws-free combination and tracks permute", {
  skip_if_not_installed("fastperm")
  skip_if_not(exists("fastperm_spa_quadratic", where = asNamespace("fastperm")),
              "fastperm lacks fastperm_spa_quadratic (needs the route-b feature)")
  ## few representations so the screen's quadratic branch is grid-feasible
  reps <- riposte_reps_default()[c("raw", "rank")]
  d <- mk_df(nb = 30L, bs = 10L, m = 5L)
  res <- riposte_test(y ~ z | block, d, statistic = "screen",
                      representations = reps, engine = "saddlepoint")
  expect_true(res$p.value > 0 && res$p.value <= 1)
  expect_true(res$combination %in% c("quadratic", "cauchy"))
  ## same calibration target as the permutation screen
  p_perm <- riposte_test(y ~ z | block, d, statistic = "screen",
                         representations = reps, engine = "permute",
                         nresample = 4999L)$p.value
  expect_lt(abs(res$p.value - p_perm), 0.08)
})

test_that("the saddlepoint screen survives a rank-deficient set instead of crashing", {
  skip_if_not_installed("fastperm")
  skip_if_not(exists("fastperm_spa_quadratic", where = asNamespace("fastperm")),
              "fastperm lacks fastperm_spa_quadratic (needs the route-b feature)")
  ## a collinear pair makes Sigma singular (the case that motivates the screen):
  ## the screen must regularize (shrink) or fall back to the Cauchy (threshold),
  ## not invert a singular metric
  reps <- list(raw = function(y) y, dup = function(y) 2 * y,
               rnk = function(y) rank(y))
  d <- mk_df(nb = 30L, bs = 10L, m = 5L)
  res <- riposte_test(y ~ z | block, d, statistic = "screen",
                      representations = reps, engine = "saddlepoint")
  expect_true(is.finite(res$p.value) && res$p.value > 0 && res$p.value <= 1)
})

test_that("the saddlepoint engine refuses max, adjustment, and too-large r", {
  skip_if_not_installed("fastperm")
  d <- mk_df()
  d$x <- rnorm(nrow(d))
  ## the max needs a joint (multivariate-max) saddlepoint not yet built
  expect_error(riposte_test(y ~ z | block, d, statistic = "max",
                            engine = "saddlepoint"), "saddlepoint")
  ## covariance adjustment is not available draws-free
  expect_error(riposte_test(y ~ z | block, d, statistic = "cauchy",
                            engine = "saddlepoint", adjust = ~ x), "saddlepoint")
  ## the quadratic saddlepoint is a tensor grid feasible only for a few
  ## representations; a high-rank set must point the user at the permute engine
  skip_if_not(exists("fastperm_spa_quadratic", where = asNamespace("fastperm")),
              "fastperm lacks fastperm_spa_quadratic (needs the route-b feature)")
  poly <- list(p1 = function(y) y,      p2 = function(y) y^2,
               p3 = function(y) y^3,    p4 = function(y) y^4,
               p5 = function(y) abs(y))
  expect_error(riposte_test(y ~ z | block, d, statistic = "quadratic",
                            representations = poly, engine = "saddlepoint"),
               "permute|representation")
})
