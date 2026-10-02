## test-exactness.R
##
## The headline property is exactness: the permutation p-value is the share of the
## randomization distribution at least as extreme as observed. The size
## simulations check this with Monte-Carlo noise; this test checks it EXACTLY. On
## a tiny design we enumerate the FULL within-block randomization distribution
## (every within-block assignment), feed it as the draws, and confirm riposte's
## p-values equal an independent brute-force computation to machine precision.
## This pins down the arithmetic exactness rests on: observed-row inclusion, the
## >= comparison, and the closed-form moments.

## enumerate_assignments() lives in helper-enumerate.R (shared with the cluster
## exactness test).

## An assignment and its within-block complement give treated-score sums T and
## -T, so |T| and the quadratic form tie EXACTLY in pairs, and floating point
## separates each pair in its last bits. The brute-force counts therefore treat
## values within 1e-9 of the observed one as ties, as riposte does.
at_least <- function(v) mean(v >= v[1] - 1e-9)

test_that("the full enumeration gives the exact permutation p-value", {
  ## 2 blocks of 6, 3 treated each -> choose(6,3)^2 = 400 assignments. (Blocks
  ## of 4 with 2 treated are too small: the 6 treated sets in a block are 3
  ## pairs +v, -v in three dimensions, so the quadratic form equals the rank of
  ## Sigma at every assignment and the exact p-value is always 1.)
  set.seed(20260623)
  block <- factor(rep(1:2, each = 6))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  expect_equal(ncol(G), 400L)

  sm <- riposte_score_matrix(y, block)
  mom <- riposte_sw_moments(sm$scores, G[, 1], block)

  ## independent brute force: the quadratic form on every enumerated assignment,
  ## using the closed-form moments, and the exact share >= observed
  Tmat <- crossprod(G, sm$scores)
  cen <- sweep(Tmat, 2, mom$mu)
  q <- rowSums((cen %*% MASS::ginv(mom$Sigma)) * cen)
  expect_gt(diff(range(q)), 1)                 # not the degenerate design
  p_exact <- at_least(q)

  ## riposte, given the full enumeration as its draws, must match exactly
  p_riposte <- riposte_quadratic(sm$scores, G[, 1], block, moments = mom, draws = G)$p.value
  expect_equal(p_riposte, p_exact)
  ## and it must be a genuine multiple of 1/400 (exact, not Monte-Carlo)
  expect_equal(p_riposte * 400, round(p_riposte * 400))
})

test_that("single-representation combinations equal the exact two-sided permutation test", {
  ## with one representation, quadratic/max/cauchy all reduce to the two-sided
  ## exact permutation test of that representation's linear statistic
  set.seed(11)
  block <- factor(rep(1:2, each = 4))
  y <- rnorm(8)
  G <- enumerate_assignments(block)
  raw <- list(raw = function(v) v)
  sm <- riposte_score_matrix(y, block, raw)

  ## independent exact two-sided p of the raw (difference-in-means) statistic
  Traw <- as.numeric(crossprod(G, sm$scores))
  p_exact <- at_least(abs(Traw))

  expect_equal(riposte_max(sm$scores, G[, 1], block, draws = G)$p.value, p_exact)
  expect_equal(riposte_quadratic(sm$scores, G[, 1], block,
                                 moments = riposte_sw_moments(sm$scores, G[, 1], block),
                                 draws = G)$p.value, p_exact)
})
