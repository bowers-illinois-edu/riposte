## test-combine.R
##
## The statistical content:
##   - The quadratic STATISTIC equals coin's quadratic energy statistic to
##     machine precision (the spec's equivalence claim). Because riposte uses the
##     closed-form moments and a Moore-Penrose inverse, it matches coin's
##     teststat = "quadratic" exactly.
##   - All three combinations return valid permutation p-values: in (0, 1], with
##     the observed assignment included so the smallest possible value is
##     1/(1 + nresample).
##   - The Cauchy combination stays finite even when a representation is sharply
##     separated --- the mid-p transform never reaches a pole (no clamp, no probit).
##   - A strong, unambiguous effect is detected by all three (a sanity check;
##     the subtle power-direction contrasts that separate the combinations are
##     statistical and live in the size/power simulations, not here).
##   - Sharing one set of draws makes the combinations reproducible.

make_design <- function(seed = 20260623, effect = 0) {
  set.seed(seed)
  block <- factor(rep(1:6, each = 10))
  z <- integer(length(block))
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(length(ix), length(ix) %/% 2L)] <- 1L
  }
  y0 <- rnorm(length(block))
  y <- y0 + effect * z                       # constant additive effect on treated
  list(y = y, z = z, block = block)
}

test_that("quadratic statistic equals coin's quadratic energy statistic", {
  skip_if_not_installed("coin")
  d <- make_design()
  sm <- riposte_score_matrix(d$y, d$block)
  mom <- riposte_sw_moments(sm$scores, d$z, d$block)
  Tobs <- crossprod(d$z, sm$scores)
  q_rip <- as.numeric((Tobs - mom$mu) %*% MASS::ginv(mom$Sigma) %*% t(Tobs - mom$mu))

  reps <- as.data.frame(lapply(riposte_reps_default(), function(f) {
    o <- numeric(length(d$y))
    for (ix in split(seq_along(d$y), d$block)) o[ix] <- f(d$y[ix]); o
  }))[, sm$kept]
  reps$trtF <- factor(d$z); reps$blockF <- d$block
  fml <- stats::as.formula(paste0("cbind(", paste(sm$kept, collapse = ", "),
                                  ") ~ trtF | blockF"))
  it <- coin::independence_test(fml, data = reps, teststat = "quadratic",
                                distribution = "asymptotic")
  expect_equal(q_rip, as.numeric(coin::statistic(it, "test")), tolerance = 1e-9)
})

test_that("all three combinations return valid permutation p-values", {
  d <- make_design()
  sm <- riposte_score_matrix(d$y, d$block)
  set.seed(1)
  for (res in list(riposte_quadratic(sm$scores, d$z, d$block, nresample = 199),
                   riposte_cauchy(sm$scores, d$z, d$block, nresample = 199),
                   riposte_max(sm$scores, d$z, d$block, nresample = 199))) {
    expect_true(res$p.value > 0 && res$p.value <= 1)
    expect_gte(res$p.value, 1 / 200)         # observed is included
    expect_true(is.finite(res$statistic))
  }
})

test_that("Cauchy stays finite under a sharply separated representation", {
  ## a perfectly separated outcome (treated all above control) pushes a
  ## representation's p-value to the extreme; mid-p keeps the transform finite
  block <- factor(rep(1:4, each = 6))
  z <- rep(c(1, 1, 1, 0, 0, 0), 4)
  y <- ifelse(z == 1, 10, 0) + rnorm(24, sd = 1e-6)
  sm <- riposte_score_matrix(y, block)
  set.seed(2)
  res <- riposte_cauchy(sm$scores, z, block, nresample = 199)
  expect_true(is.finite(res$statistic))
  expect_true(res$p.value > 0 && res$p.value <= 1)
})

test_that("a strong constant effect is detected by all three", {
  d <- make_design(effect = 3)               # large shift: unambiguous signal
  sm <- riposte_score_matrix(d$y, d$block)
  set.seed(3)
  expect_lt(riposte_quadratic(sm$scores, d$z, d$block, nresample = 499)$p.value, 0.05)
  expect_lt(riposte_cauchy(sm$scores, d$z, d$block, nresample = 499)$p.value, 0.05)
  expect_lt(riposte_max(sm$scores, d$z, d$block, nresample = 499)$p.value, 0.05)
})

test_that("max identifies a location-sensitive representation for a location shift", {
  ## A symmetric location shift loads on the monotone location representations
  ## (raw, rank, huber) and is nearly invisible to the distance representations
  ## (treated and control are equally far apart, so their centred distance scores
  ## roughly cancel). This is exactly why combining several representations helps:
  ## no single representation sees every kind of effect.
  d <- make_design(effect = 3)
  sm <- riposte_score_matrix(d$y, d$block)
  set.seed(4)
  res <- riposte_max(sm$scores, d$z, d$block, nresample = 199)
  expect_true(res$which %in% c("raw", "rank", "huber"))
  expect_false(res$which %in% c("mean_dist", "mean_rank_dist", "max_dist"))
})

test_that("sharing one set of draws makes combinations reproducible", {
  d <- make_design(effect = 1)
  sm <- riposte_score_matrix(d$y, d$block)
  set.seed(5); G <- riposte_block_draws(d$z, d$block, nresample = 299)
  a <- riposte_quadratic(sm$scores, d$z, d$block, draws = G)
  b <- riposte_quadratic(sm$scores, d$z, d$block, draws = G)
  expect_identical(a$p.value, b$p.value)
  ## and a seeded from-scratch run reproduces a shared-draws run
  set.seed(6); c1 <- riposte_cauchy(sm$scores, d$z, d$block, nresample = 299)
  set.seed(6); G2 <- riposte_block_draws(d$z, d$block, nresample = 299)
  c2 <- riposte_cauchy(sm$scores, d$z, d$block, draws = G2)
  expect_equal(c1$p.value, c2$p.value)
})
