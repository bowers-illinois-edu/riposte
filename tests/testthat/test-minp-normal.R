## test-minp-normal.R
##
## THE POINT. Jake's task of 2026-11-06 (and the second of the three ways of
## computing the min-p in the memo of decisions of 12 November): draw the
## standardized statistics from their joint normal approximation instead of
## redrawing assignments, compute every p-value of the min-p from each draw, and
## refer the observed smallest to the smallest at the draws. The six
## standardized sums are N(0, R) with R their null correlation matrix; the
## quadratic combination is a function of the same six; and each block statistic
## is a sum of squared per-block standardized sums, which are jointly normal with
## the six, with covariances from the same within-block moments. The work per
## draw does not grow with the number of units. Under the normal model each
## p-value from the normal or chi-square formula is exactly uniform, so the
## formula p-values are the ones riposte uses here. riposte_test(engine =
## "asymptotic") now runs this min-p, with nresample normal draws and `seed`.
## These tests check:
##   1. the joint covariance of the treated sums and the per-block sums equals
##      their covariance over every assignment the design allows;
##   2. the observed p-values are the formula p-values the permutation engine
##      reports under pvalue = "asymptotic";
##   3. with the six single p-values only, the large-sample min-p is the
##      large-sample max-T, up to Monte Carlo error;
##   4. in a large experiment it is close to the permutation min-p;
##   5. it is reproducible with `seed`, prints that it is approximate, and
##      leaves the quadratic combination and block statistics out of a
##      one-sided test.

normal_design <- function(seed, B = 10L, nb = 40L, effect = 0.2) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  block <- factor(rep(seq_len(B), each = nb))
  z <- integer(B * nb)
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L
  }
  data.frame(y = rnorm(B * nb, sd = ifelse(z == 1, 1 + effect, 1)), z = z, b = block)
}

test_that("the joint covariance of the treated sums and the per-block sums is their covariance over every assignment", {
  set.seed(20261118)
  block <- factor(rep(1:3, each = 4))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  sm <- riposte_score_matrix(y, block)
  rank_b <- within_block_rank(y, block)
  ## the treated sums of the six scores and, for the rank, each block's treated
  ## sum of departures from its block mean
  blocks_dev <- sapply(levels(block), function(b) {
    a <- ifelse(block == b, rank_b - stats::ave(rank_b, block), 0)
    a
  })
  all_sums <- crossprod(G, cbind(sm$scores, blocks_dev))
  expected <- enum_moments(all_sums)$Sigma
  got <- riposte:::riposte_joint_null_cov(sm$scores, G[, 1], block, cbind(rank = rank_b))
  expect_equal(unname(got$Sigma), unname(expected), tolerance = 1e-10)
})

test_that("the large-sample min-p reports the formula p-values the permutation engine reports", {
  d <- normal_design(1)
  ctl <- riposte_minp_control(block = "rank")
  asy <- riposte_test(y ~ z | b, d, engine = "asymptotic", nresample = 999, seed = 1,
                      minp = ctl)
  per <- riposte_test(y ~ z | b, d, nresample = 199, seed = 1,
                      minp = riposte_minp_control(pvalue = "asymptotic", block = "rank"))
  expect_equal(asy$component_p, per$component_p, tolerance = 1e-8)
  expect_equal(asy$requested, "minp")
  expect_equal(length(asy$component_p), 8L)
})

test_that("with the six single p-values only, the large-sample min-p is the large-sample max-T", {
  skip_if_not_installed("mvtnorm")
  d <- normal_design(2)
  minp <- riposte_test(y ~ z | b, d, statistic = "minp", engine = "asymptotic",
                       nresample = 19999, seed = 2,
                       minp = riposte_minp_control(quadratic = FALSE))
  maxt <- riposte_test(y ~ z | b, d, statistic = "max", engine = "asymptotic")
  se <- sqrt(maxt$p.value * (1 - maxt$p.value) / 20000)
  expect_lt(abs(minp$p.value - maxt$p.value), 4 * se + 1e-4)
})

test_that("in a large experiment the large-sample min-p is close to the permutation min-p", {
  skip_on_cran()
  d <- normal_design(3, B = 20L, nb = 50L)
  asy <- riposte_test(y ~ z | b, d, engine = "asymptotic", nresample = 9999, seed = 3)
  per <- riposte_test(y ~ z | b, d, nresample = 4999, seed = 3,
                      minp = riposte_minp_control(pvalue = "asymptotic"))
  expect_lt(abs(asy$p.value - per$p.value), 0.03)
})

test_that("the large-sample min-p is reproducible, says it is approximate, and is one-sided without the squared statistics", {
  d <- normal_design(4)
  a <- riposte_test(y ~ z | b, d, engine = "asymptotic", nresample = 499, seed = 4)
  b <- riposte_test(y ~ z | b, d, engine = "asymptotic", nresample = 499, seed = 4)
  expect_equal(a$p.value, b$p.value)
  expect_equal(a$nresample, 499L)
  expect_output(print(a), "approximate")
  expect_output(print(a), "normal")
  g <- riposte_test(y ~ z | b, d, engine = "asymptotic", alternative = "greater",
                    nresample = 499, seed = 4)
  expect_equal(length(g$component_p), 6L)
  expect_error(riposte_test(y ~ z | b, d, engine = "asymptotic", alternative = "greater",
                            minp = riposte_minp_control(block = "rank")),
               "both directions")
})
