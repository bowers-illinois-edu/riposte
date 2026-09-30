## test-moments.R
##
## The statistical content: riposte computes the EXACT permutation mean and
## covariance of the representation linear statistics in closed form
## (Strasser-Weber 1999), without re-randomizing. Three independent checks:
##   1. The covariance matches coin's covariance() to machine precision (coin
##      implements the same linear-statistic framework). This is what makes
##      riposte's quadratic match coin's quadratic energy column.
##   2. The closed-form mean and covariance match a Monte-Carlo estimate from the
##      actual within-block re-randomizations, within sampling error. This
##      validates the formula without reference to coin.
##   3. Structure: Sigma is symmetric and positive semidefinite; a block with no
##      randomization (all treated or all control) contributes nothing; the mean
##      of a within-block-centred statistic is zero.

## a small unbalanced block design reused across tests
make_design <- function(seed = 20260623) {
  set.seed(seed)
  block <- factor(rep(1:4, times = c(6, 8, 4, 10)))
  y <- rnorm(length(block))
  z <- integer(length(block))
  for (b in levels(block)) {
    ix <- which(block == b)
    z[ix][sample.int(length(ix), length(ix) %/% 2L)] <- 1L
  }
  list(y = y, z = z, block = block)
}

test_that("closed-form covariance matches coin to machine precision (univariate)", {
  skip_if_not_installed("coin")
  d <- make_design()
  ## raw representation only; coin computes the same linear statistic
  scores <- matrix(d$y, ncol = 1, dimnames = list(NULL, "raw"))
  m <- riposte_sw_moments(scores, d$z, d$block)
  dat <- data.frame(Y = d$y, trtF = factor(d$z), blockF = d$block)
  it <- coin::independence_test(Y ~ trtF | blockF, data = dat,
                                distribution = "asymptotic")
  expect_equal(as.numeric(m$Sigma), as.numeric(coin::covariance(it)),
               tolerance = 1e-10)
})

test_that("closed-form covariance matches coin for all six representations", {
  skip_if_not_installed("coin")
  d <- make_design()
  sm <- riposte_score_matrix(d$y, d$block)             # the six reps, centred
  m <- riposte_sw_moments(sm$scores, d$z, d$block)
  ## build coin on the SAME representation values (raw, uncentred); the
  ## covariance is centring-invariant, so coin's matrix must equal ours
  reps <- as.data.frame(lapply(riposte_reps_default(), function(f) {
    out <- numeric(length(d$y))
    for (ix in split(seq_along(d$y), d$block)) out[ix] <- f(d$y[ix])
    out
  }))
  reps <- reps[, sm$kept, drop = FALSE]
  reps$trtF <- factor(d$z); reps$blockF <- d$block
  fml <- stats::as.formula(paste0("cbind(", paste(sm$kept, collapse = ", "),
                                  ") ~ trtF | blockF"))
  it <- coin::independence_test(fml, data = reps, teststat = "quadratic",
                                distribution = "asymptotic")
  expect_equal(unname(m$Sigma), unname(drop(coin::covariance(it))),
               tolerance = 1e-8)
})

test_that("closed-form mean and covariance match a Monte-Carlo estimate", {
  d <- make_design()
  sm <- riposte_score_matrix(d$y, d$block)
  m <- riposte_sw_moments(sm$scores, d$z, d$block)
  ## the randomization distribution of T = sum_treated (centred score), drawn the
  ## way the study assigns treatment; drop column 1 (the observed assignment)
  set.seed(7)
  G <- riposte_block_draws(d$z, d$block, nresample = 20000)[, -1, drop = FALSE]
  Tdraws <- crossprod(G, sm$scores)                    # draws x representations
  emp_mu <- colMeans(Tdraws)
  emp_Sigma <- stats::cov(Tdraws)
  ## centred scores => permutation mean zero; MC mean is near zero
  expect_true(max(abs(emp_mu)) < 0.15)
  expect_equal(m$mu, rep(0, length(m$mu)), ignore_attr = TRUE, tolerance = 1e-8)
  ## covariance matches within Monte-Carlo error (relative)
  expect_equal(emp_Sigma, m$Sigma, ignore_attr = TRUE, tolerance = 0.06)
})

test_that("Sigma is symmetric and positive semidefinite", {
  d <- make_design()
  sm <- riposte_score_matrix(d$y, d$block)
  m <- riposte_sw_moments(sm$scores, d$z, d$block)
  expect_equal(m$Sigma, t(m$Sigma))
  ev <- eigen(m$Sigma, symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(ev) > -1e-8)
})

test_that("a block with no randomization contributes nothing to Sigma", {
  ## block 2 is all treated: m = n, so w = m(n-m)/(n(n-1)) = 0, no contribution
  y <- c(rnorm(6), rnorm(6))
  block <- factor(rep(1:2, each = 6))
  z <- c(rep(c(1, 0), 3), rep(1, 6))                   # block 2 all treated
  scores <- matrix(y, ncol = 1, dimnames = list(NULL, "raw"))
  m_full <- riposte_sw_moments(scores, z, block)
  ## same as keeping only block 1
  m_b1 <- riposte_sw_moments(scores[block == 1, , drop = FALSE],
                             z[block == 1], droplevels(block[block == 1]))
  expect_equal(as.numeric(m_full$Sigma), as.numeric(m_b1$Sigma))
})

test_that("the mean of a within-block-centred statistic is zero", {
  d <- make_design()
  sm <- riposte_score_matrix(d$y, d$block)             # centred scores
  m <- riposte_sw_moments(sm$scores, d$z, d$block)
  expect_true(all(abs(m$mu) < 1e-10))
})
