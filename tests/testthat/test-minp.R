## test-minp.R
##
## THE POINT. Single-step min-p is the fifth combination. On every assignment
## it takes each representation's mid-p value and keeps the smallest; its
## p-value is the share of assignments whose smallest mid-p value is at most
## the observed one. These tests pin down four things the combination rests on:
##   1. exactness: on a full enumeration the p-value equals a brute-force count;
##   2. the identity with the max: when every representation's statistic has
##      the same null distribution, min-p and the max give the same p-value,
##      and with one representation both are that representation's own test;
##   3. where they part: a representation whose statistic takes few values
##      makes min-p and the max disagree;
##   4. the Bonferroni bracket: min-p's p-value is never below the smallest
##      observed mid-p value and never above about K times it.
## Then the level under the sharp null, including a two-valued representation,
## and the riposte_test() interface.

## enumerate_assignments() lives in helper-enumerate.R.

## Brute-force min-p on a full enumeration: Tmat has one row per assignment,
## row 1 observed. Mid-p values come from the same function the Cauchy
## combination uses, so the test checks the combination, not the mid-p
## arithmetic (test-pvalues.R does that).
brute_minp <- function(Tmat, alternative = "two.sided") {
  P <- apply(Tmat, 2, riposte_midp, alternative = alternative)
  smallest <- apply(P, 1, min)
  ## the smallest mid-p values tie exactly in pairs (an assignment and its
  ## within-block complement give sums T and -T); count values within 1e-9
  ## as ties, as riposte does
  list(p = mean(smallest <= smallest[1] + 1e-9), smallest = smallest, P = P)
}

## the max's statistic on the same enumeration, for the identity tests
brute_max_stat <- function(Tmat, mom, alternative = "two.sided") {
  Z <- sweep(sweep(Tmat, 2, mom$mu), 2, sqrt(diag(mom$Sigma)), "/")
  Zdir <- switch(alternative, two.sided = abs(Z), greater = Z, less = -Z)
  apply(Zdir, 1, max)
}

## the top indicator: 1 for the largest outcome in each block, 0 otherwise. Its
## treated sum counts the blocks whose top unit is treated, so with two blocks
## it takes only the values 0, 1, 2. This is the few-valued representation that
## separates min-p from the max.
top_in_block <- function(v) as.numeric(v == max(v))

## ---- exactness ---------------------------------------------------------------

test_that("min-p equals a brute-force count on the full enumeration, two-sided", {
  ## 2 blocks of 6, 3 treated each -> 400 assignments
  set.seed(20261005)
  block <- factor(rep(1:2, each = 6))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  sm <- riposte_score_matrix(y, block)
  Tmat <- crossprod(G, sm$scores)
  colnames(Tmat) <- colnames(sm$scores)
  bf <- brute_minp(Tmat)

  res <- riposte:::riposte_minp(sm$scores, G[, 1], block, draws = G)
  expect_equal(res$p.value, bf$p)
  ## exact, not Monte-Carlo: a multiple of 1/400
  expect_equal(res$p.value * 400, round(res$p.value * 400))
  ## the statistic is the observed smallest mid-p value, and the function names
  ## the representation that gave it
  expect_equal(res$statistic, bf$smallest[1])
  expect_equal(res$which, colnames(sm$scores)[which.min(bf$P[1, ])])
  ## the six observed mid-p values come back, as the hybrid's do
  expect_equal(res$component_p, stats::setNames(bf$P[1, ], colnames(sm$scores)))
})

test_that("one-sided min-p equals a brute-force count on the full enumeration", {
  set.seed(20261002)
  block <- factor(rep(1:2, each = 6))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  reps <- riposte_poly_reps(c(2, 5), tail = "lower")
  sm <- riposte_score_matrix(y, block, reps)
  Tmat <- crossprod(G, sm$scores)
  for (alt in c("greater", "less")) {
    bf <- brute_minp(Tmat, alt)
    res <- riposte:::riposte_minp(sm$scores, G[, 1], block, draws = G, alternative = alt)
    expect_equal(res$p.value, bf$p)
    expect_equal(res$p.value * 400, round(res$p.value * 400))
  }
})

## ---- the identity with the max -------------------------------------------------

test_that("with one representation, min-p is the exact two-sided permutation test", {
  ## the mid-p value is a decreasing function of |T|, so ordering assignments
  ## by their smallest (only) mid-p value is ordering them by |T|
  set.seed(11)
  block <- factor(rep(1:2, each = 4))
  y <- rnorm(8)
  G <- enumerate_assignments(block)
  raw <- list(raw = function(v) v)
  sm <- riposte_score_matrix(y, block, raw)
  Traw <- as.numeric(crossprod(G, sm$scores))
  p_exact <- mean(abs(Traw) >= abs(Traw[1]) - 1e-9)
  expect_equal(riposte:::riposte_minp(sm$scores, G[, 1], block, draws = G)$p.value, p_exact)
})

test_that("when every statistic has the same null distribution, min-p equals the max", {
  ## Two representations: the outcome, and the same outcome values shuffled
  ## among the units within each block. Their treated sums draw from the same
  ## values in each block, so they have the same null distribution, but they
  ## take different values on each assignment. The smallest mid-p value is then
  ## a decreasing function of the largest |z|, so the two combinations order the
  ## 400 assignments the same way and give the same p-value, whichever
  ## assignment is the observed one.
  set.seed(20261005)
  block <- factor(rep(1:2, each = 6))
  y <- rnorm(12)
  y_shuf <- y
  for (b in levels(block)) { ix <- which(block == b); y_shuf[ix] <- sample(y[ix]) }
  G <- enumerate_assignments(block)
  ## centre each column within blocks the way riposte does
  sm <- riposte_score_matrix(y, block, list(a = function(v) v))
  sm2 <- riposte_score_matrix(y_shuf, block, list(b = function(v) v))
  S <- cbind(a = sm$scores[, 1], b = sm2$scores[, 1])
  mom <- riposte_sw_moments(S, G[, 1], block)
  ## same mean and standard deviation, as the identity needs
  expect_equal(mom$mu[["a"]], mom$mu[["b"]])
  expect_equal(mom$Sigma["a", "a"], mom$Sigma["b", "b"])

  Tmat <- crossprod(G, S)
  colnames(Tmat) <- c("a", "b")
  smallest <- brute_minp(Tmat)$smallest
  largest <- brute_max_stat(Tmat, mom)
  ## for every assignment as the observed one, the two p-values agree
  p_minp <- vapply(seq_len(ncol(G)), function(j) mean(smallest <= smallest[j] + 1e-9), 0)
  p_max  <- vapply(seq_len(ncol(G)), function(j) mean(largest >= largest[j] - 1e-9), 0)
  expect_equal(p_minp, p_max)
  ## and riposte's two functions agree on the observed assignment
  expect_equal(riposte:::riposte_minp(S, G[, 1], block, moments = mom, draws = G)$p.value,
               riposte:::riposte_max(S, G[, 1], block, moments = mom, draws = G)$p.value)
})

test_that("a few-valued representation makes min-p and the max disagree", {
  ## the top indicator's treated sum takes only three values, so a standardized
  ## value that is rare for the raw outcome can be common for it; the two
  ## combinations then rank some assignments differently
  set.seed(20261005)
  block <- factor(rep(1:2, each = 6))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  sm <- riposte_score_matrix(y, block, list(raw = function(v) v, top = top_in_block))
  mom <- riposte_sw_moments(sm$scores, G[, 1], block)
  Tmat <- crossprod(G, sm$scores)
  colnames(Tmat) <- colnames(sm$scores)
  smallest <- brute_minp(Tmat)$smallest
  largest <- brute_max_stat(Tmat, mom)
  p_minp <- vapply(seq_len(ncol(G)), function(j) mean(smallest <= smallest[j] + 1e-9), 0)
  p_max  <- vapply(seq_len(ncol(G)), function(j) mean(largest >= largest[j] - 1e-9), 0)
  differ <- which(abs(p_minp - p_max) > 1e-9)
  expect_gt(length(differ), 0)
  ## riposte reproduces the disagreement with that assignment observed
  j <- differ[1]
  Gj <- cbind(G[, j], G[, -j])
  expect_false(isTRUE(all.equal(
    riposte:::riposte_minp(sm$scores, G[, j], block, moments = mom, draws = Gj)$p.value,
    riposte:::riposte_max(sm$scores, G[, j], block, moments = mom, draws = Gj)$p.value)))
})

## ---- the Bonferroni bracket ----------------------------------------------------

test_that("min-p's p-value is never below the smallest observed mid-p value", {
  ## The share of assignments whose smallest mid-p value is at most the observed
  ## smallest, c, is at least the share whose mid-p value for that same
  ## representation is at most c, and that share is at least c because the
  ## mid-p value counts ties as half while the share counts them in full.
  set.seed(20261005)
  block <- factor(rep(1:2, each = 6))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  sm <- riposte_score_matrix(y, block)
  res <- riposte:::riposte_minp(sm$scores, G[, 1], block, draws = G)
  expect_gte(res$p.value, res$statistic - 1e-12)
})

test_that("min-p's p-value is never above the Bonferroni-adjusted p-value, up to tie slack", {
  ## One-sided, with three continuous representations, so no two assignments
  ## give the same sum: each representation's mid-p values are then the n
  ## numbers (j - 0.5) / n. The share of assignments whose mid-p value for one
  ## representation is at most c is at most c + 0.5 / n, and the share whose
  ## smallest mid-p value is at most c is at most the sum of those K shares.
  set.seed(20261005)
  block <- factor(rep(1:2, each = 6))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  reps <- list(raw = function(v) v, cube = function(v) v^3, expo = function(v) exp(v))
  sm <- riposte_score_matrix(y, block, reps)
  K <- ncol(sm$scores); n <- ncol(G)
  for (alt in c("greater", "less")) {
    res <- riposte:::riposte_minp(sm$scores, G[, 1], block, draws = G, alternative = alt)
    expect_lte(res$p.value, K * res$statistic + K * 0.5 / n + 1e-12)
    expect_gte(res$p.value, res$statistic - 1e-12)
  }
})

## ---- level under the sharp null ------------------------------------------------

test_that("under the sharp null min-p holds the 0.05 level, with a two-valued representation", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261005)
  B <- 6L; nb <- 10L; nsims <- 300L; nres <- 199L; alpha <- 0.05
  block <- factor(rep(seq_len(B), each = nb))
  reps <- c(riposte_reps_default(), list(top = top_in_block))
  rej <- matrix(0, nsims, 2, dimnames = list(NULL, c("two.sided", "greater")))
  for (i in seq_len(nsims)) {
    y <- rnorm(B * nb)                 # sharp null: no effect anywhere
    z <- integer(B * nb)
    for (b in levels(block)) {
      ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L
    }
    sm <- riposte_score_matrix(y, block, reps)
    G <- riposte_block_draws(z, block, nres)
    for (alt in colnames(rej))
      rej[i, alt] <- riposte:::riposte_minp(sm$scores, z, block, draws = G,
                                             alternative = alt)$p.value <= alpha
  }
  ## SE ~ sqrt(0.05 * 0.95 / 300) = 0.0126; the band is about +/- 3 SE, as in
  ## test-size.R
  for (nm in colnames(rej)) {
    expect_lt(mean(rej[, nm]), 0.09)
    expect_gt(mean(rej[, nm]), 0.015)
  }
})

## ---- the riposte_test() interface ----------------------------------------------

minp_design <- function(seed, B = 4L, nb = 8L) {
  set.seed(seed)
  block <- factor(rep(seq_len(B), each = nb))
  z <- integer(B * nb)
  for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L }
  data.frame(y = rnorm(B * nb) + 0.8 * z, z = z, b = block)
}

test_that("riposte_test(statistic = 'minp') runs and reports what the max reports", {
  dat <- minp_design(1)
  res <- riposte_test(y ~ z | b, dat, statistic = "minp", nresample = 199, seed = 1)
  expect_s3_class(res, "riposte_test")
  expect_equal(res$combination, "minp")
  expect_equal(res$requested, "minp")
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_equal(length(res$component_p), 6L)
  expect_true(all(res$component_p > 0 & res$component_p < 1))
  expect_equal(res$statistic, min(res$component_p))
  expect_output(print(res), "minp")
})

test_that("min-p accepts a one-sided alternative", {
  dat <- minp_design(2)
  res <- riposte_test(y ~ z | b, dat, statistic = "minp", alternative = "greater",
                      nresample = 199, seed = 2)
  expect_equal(res$combination, "minp")
  expect_equal(res$alternative, "greater")
  ## the effect raises treated outcomes, so "greater" finds it and "less" does not
  res_less <- riposte_test(y ~ z | b, dat, statistic = "minp", alternative = "less",
                           nresample = 199, seed = 2)
  expect_lt(res$p.value, res_less$p.value)
})

test_that("min-p runs with covariance adjustment", {
  dat <- minp_design(3)
  dat$x <- dat$y + rnorm(nrow(dat))
  res <- riposte_test(y ~ z | b, dat, statistic = "minp", adjust = ~ x,
                      nresample = 99, seed = 3)
  expect_equal(res$combination, "minp")
  expect_true(res$adjusted)
  expect_true(res$p.value >= 0 && res$p.value <= 1)
})

test_that("the large-sample and saddlepoint engines refuse min-p", {
  dat <- minp_design(4)
  expect_error(riposte_test(y ~ z | b, dat, statistic = "minp", engine = "asymptotic"),
               "minp|min-p")
  expect_error(riposte_test(y ~ z | b, dat, statistic = "minp", engine = "saddlepoint"),
               "minp|min-p")
})

test_that("riposte_test() requires the full name of the statistic", {
  ## "max" and "minp" share a first letter, so riposte_test() accepts only full
  ## names; an abbreviation stops with an error that lists the names
  dat <- minp_design(5)
  for (abbrev in c("m", "mi", "ma", "quad"))
    expect_error(riposte_test(y ~ z | b, dat, statistic = abbrev, nresample = 99),
                 "spelled in full")
  expect_equal(riposte_test(y ~ z | b, dat, statistic = "max", nresample = 99,
                            seed = 5)$combination, "max")
  expect_equal(riposte_test(y ~ z | b, dat, statistic = "minp", nresample = 99,
                            seed = 5)$combination, "minp")
  ## the default is unchanged
  expect_equal(riposte_test(y ~ z | b, dat, nresample = 99, seed = 5)$requested, "screen")
})
