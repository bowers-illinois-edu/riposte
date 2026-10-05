## test-one-sided.R
##
## The statistical content:
##   - The bottom-weighted polynomial score of a unit at ascending within-block
##     rank k of n_b is -((n_b + 1 - k) / (n_b + 1))^(zeta - 1). It is defined
##     with a minus sign so that it RISES with the outcome, like the raw outcome,
##     the rank, and the top-weighted score. Then alternative = "less" means
##     "treated outcomes are lower" for every one of these representations at
##     once, which is how a user asks about harm.
##   - A one-sided test counts only evidence in the named direction. With one
##     representation, its p-value is the exact one-sided permutation p-value of
##     the treated-score sum; with several, the max combination takes the largest
##     standardized sum in that direction, and the Cauchy combination averages
##     the Cauchy transforms of one-sided mid-p values.
##   - Every one-sided test is a test of the sharp null of no effect, so it holds
##     its level under that null.
##   - With rank-based representations and no ties, the one-sided test is also
##     valid for the hypothesis that no unit was harmed (Caughey, Dafoe, Li, and
##     Miratrix 2023). The core of that argument is deterministic: if every
##     treated unit's outcome is at least its no-effect value, then at EVERY
##     assignment the "less" p-value computed from the observed outcomes is at
##     least the one computed from the no-effect outcomes. We check that on a
##     full enumeration.
##   - Against harm to a few treated units, the one-sided test with
##     bottom-weighted scores rejects more often than the two-sided difference
##     in means.
##   - The quadratic and the screen have no one-sided form, so they refuse a
##     one-sided alternative, and a one-sided call that names no statistic uses
##     the max. The saddlepoint and asymptotic engines compute the one-sided
##     Cauchy from one-sided marginal p-values.

## a small block design with a fixed seed; tau is added to the treated outcomes
one_sided_design <- function(seed, B = 6L, nb = 10L, tau = 0) {
  set.seed(seed)
  block <- factor(rep(seq_len(B), each = nb))
  z <- integer(B * nb)
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L
  }
  y0 <- rnorm(B * nb)
  list(y0 = y0, y = y0 + tau * z, z = z, block = block)
}

## ---- the bottom-weighted score ---------------------------------------------

test_that("the bottom-weighted score is minus ((n + 1 - k) / (n + 1))^(zeta - 1)", {
  y <- c(3.1, -0.4, 2.2, 2.2, 7.5, 0.0)
  k <- rank(y, ties.method = "average")               # the two 2.2s share 3.5
  n <- length(y)
  expect_equal(riposte_poly_scores(y, 5, tail = "lower"),
               -((n + 1 - k) / (n + 1))^4)
  ## without ties, ranking from the top is ranking -y from the bottom
  y2 <- c(3.1, -0.4, 2.2, 7.5, 0.0)
  expect_equal(riposte_poly_scores(y2, 5, tail = "lower"),
               -riposte_poly_scores(-y2, 5))
  ## the default tail is unchanged
  expect_equal(riposte_poly_scores(y, 5), (k / (n + 1))^4)
  expect_equal(riposte_poly_scores(y, 5, tail = "upper"), (k / (n + 1))^4)
  expect_error(riposte_poly_scores(y, 5, tail = "middle"))
})

test_that("both tails rise with the outcome, and the lower one weights the bottom", {
  y <- c(5, 1, 9, 3, 7, 2, 8)
  ord <- order(y)
  for (tail in c("upper", "lower")) {
    s <- riposte_poly_scores(y, 14, tail = tail)
    expect_true(all(diff(s[ord]) > 0))
  }
  ## the lower tail spends its range at the bottom: the step from the lowest to
  ## the second-lowest unit is larger than the step from the second-highest to
  ## the highest; the upper tail is the reverse
  low <- riposte_poly_scores(y, 14, tail = "lower")[ord]
  up  <- riposte_poly_scores(y, 14, tail = "upper")[ord]
  expect_gt(low[2] - low[1], low[7] - low[6])
  expect_gt(up[7] - up[6], up[2] - up[1])
})

test_that("riposte_poly_reps names the lower tail apart from the upper tail", {
  up  <- riposte_poly_reps(c(5, 14))
  low <- riposte_poly_reps(c(5, 14), tail = "lower")
  expect_named(up, c("poly5", "poly14"))
  expect_named(low, c("polylow5", "polylow14"))
  y <- c(4, 1, 3, 2)
  expect_equal(low$polylow14(y), riposte_poly_scores(y, 14, tail = "lower"))
  ## both sets together have unique names, so they can be combined
  expect_silent(riposte:::riposte_validate_reps(c(up, low)))
})

## ---- one-sided mid-p values ---------------------------------------------------

test_that("one-sided mid-p values count only the named direction", {
  v <- c(2, -3, 0.5, 1, -1)
  n <- length(v)
  ## "greater": the share at or above each value, half the tie mass at it
  expect_equal(riposte_midp(v, "greater"), (n - rank(v) + 0.5) / n)
  expect_equal(riposte_midp(v, "less"), (n - rank(-v) + 0.5) / n)
  ## the largest value gets the smallest "greater" p and the largest "less" p
  expect_equal(which.min(riposte_midp(v, "greater")), which.max(v))
  expect_equal(which.min(riposte_midp(v, "less")), which.min(v))
  ## with no ties the two one-sided mid-p values add to one
  expect_equal(riposte_midp(v, "greater") + riposte_midp(v, "less"), rep(1, n))
  ## the default is unchanged
  expect_equal(riposte_midp(v), riposte_midp(v, "two.sided"))
})

## ---- exactness on a full enumeration ------------------------------------------

test_that("with one representation, the one-sided tests are the exact one-sided permutation tests", {
  ## 2 blocks of 4, 2 treated each: 36 assignments
  set.seed(11)
  block <- factor(rep(1:2, each = 4))
  y <- rnorm(8)
  G <- enumerate_assignments(block)
  raw <- list(raw = function(v) v)
  sm <- riposte_score_matrix(y, block, raw)
  Traw <- as.numeric(crossprod(G, sm$scores))
  mom <- riposte_sw_moments(sm$scores, G[, 1], block)

  expect_equal(riposte:::riposte_max(sm$scores, G[, 1], block, moments = mom,
                                     draws = G, alternative = "greater")$p.value,
               mean(Traw >= Traw[1] - 1e-9))
  expect_equal(riposte:::riposte_max(sm$scores, G[, 1], block, moments = mom,
                                     draws = G, alternative = "less")$p.value,
               mean(Traw <= Traw[1] + 1e-9))
})

test_that("the one-sided max equals a brute-force count on the full enumeration", {
  set.seed(20261002)
  block <- factor(rep(1:2, each = 6))                 # 400 assignments
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  reps <- riposte_poly_reps(c(2, 5), tail = "lower")
  sm <- riposte_score_matrix(y, block, reps)
  mom <- riposte_sw_moments(sm$scores, G[, 1], block)

  ## independent brute force: standardize every sum by its exact standard
  ## deviation, take the most negative one at each assignment, count
  std <- sweep(crossprod(G, sm$scores), 2, sqrt(diag(mom$Sigma)), "/")
  smallest <- apply(std, 1, min)
  ## rank sums tie exactly across assignments; floating point separates the
  ## ties in their last bits, so count values within 1e-9 as ties
  p_exact <- mean(smallest <= smallest[1] + 1e-9)

  res <- riposte:::riposte_max(sm$scores, G[, 1], block, moments = mom,
                               draws = G, alternative = "less")
  expect_equal(res$p.value, p_exact)
  expect_equal(res$p.value * 400, round(res$p.value * 400))
})

test_that("if no treated unit is harmed, the 'less' p-value never falls below its no-effect value", {
  ## The deterministic heart of the bounded-null argument. y0 holds the
  ## outcomes under control. At each assignment z, some treated units are
  ## helped (tau >= 0) and none harmed. With rank scores and no ties, the
  ## p-value computed from the observed outcomes y0 + tau * z must be at least
  ## the p-value computed from y0, at EVERY one of the 400 assignments.
  set.seed(7)
  block <- factor(rep(1:2, each = 6))
  y0 <- rnorm(12)
  tau <- c(0, 2.5, 0, 0.7, 0, 0, 1.1, 0, 0, 3, 0, 0.2)
  G <- enumerate_assignments(block)
  reps <- c(list(rank = riposte_rank), riposte_poly_reps(c(5, 14), tail = "lower"))

  p_less <- function(y, z) {
    ## reorder the enumeration so the assignment being tested comes first
    own <- which(colSums(G != z) == 0)
    draws <- cbind(G[, own], G[, -own])
    sm <- riposte_score_matrix(y, block, reps)
    mom <- riposte_sw_moments(sm$scores, z, block)
    c(max = riposte:::riposte_max(sm$scores, z, block, moments = mom,
                                  draws = draws, alternative = "less")$p.value,
      cauchy = riposte:::riposte_cauchy(sm$scores, z, block, draws = draws,
                                        alternative = "less")$p.value)
  }
  for (j in seq_len(ncol(G))) {
    z <- G[, j]
    p_obs <- p_less(y0 + tau * z, z)
    p_null <- p_less(y0, z)
    expect_true(all(p_obs >= p_null - 1e-12))
  }
})

## ---- the riposte_test interface -----------------------------------------------

test_that("a one-sided call with no statistic uses the max", {
  d <- one_sided_design(1)
  dat <- data.frame(y = d$y, z = d$z, b = d$block)
  res <- riposte_test(y ~ z | b, dat, alternative = "less",
                      representations = riposte_poly_reps(c(2, 5), tail = "lower"),
                      nresample = 99, seed = 1)
  expect_equal(res$combination, "max")
  expect_equal(res$alternative, "less")
  expect_output(print(res), "less")
  ## the two-sided default is unchanged
  res2 <- riposte_test(y ~ z | b, dat, nresample = 99, seed = 1)
  expect_equal(res2$alternative, "two.sided")
  expect_equal(res2$requested, "screen")
})

test_that("the quadratic and the screen refuse a one-sided alternative", {
  d <- one_sided_design(2)
  dat <- data.frame(y = d$y, z = d$z, b = d$block)
  expect_error(riposte_test(y ~ z | b, dat, alternative = "less",
                            statistic = "quadratic", nresample = 99), "one-sided")
  expect_error(riposte_test(y ~ z | b, dat, alternative = "less",
                            statistic = "screen", nresample = 99), "one-sided")
})

test_that("the saddlepoint Cauchy uses one-sided saddlepoint marginals", {
  skip_if_not_installed("fastperm")
  d <- one_sided_design(6, B = 20L, nb = 10L)
  dat <- data.frame(y = d$y, z = d$z, b = d$block)
  low <- riposte_poly_reps(c(2, 5), tail = "lower")
  res <- riposte_test(y ~ z | b, dat, representations = low, statistic = "cauchy",
                      engine = "saddlepoint", alternative = "less")
  ## each marginal is fastperm's one-sided "less" tail of the treated sum
  sm <- riposte_score_matrix(d$y, d$block, low)
  p_each <- vapply(seq_len(ncol(sm$scores)), function(j)
    fastperm::fastperm_spa_linear(sm$scores[, j], d$z, d$block,
                                  alternative = "less")$p.value, numeric(1))
  expect_equal(unname(res$component_p), p_each)
  expect_equal(res$p.value, 0.5 - atan(mean(riposte_acat_term(p_each))) / pi)
  ## and it tracks the permutation one-sided Cauchy
  p_perm <- riposte_test(y ~ z | b, dat, representations = low, statistic = "cauchy",
                         alternative = "less", nresample = 4999, seed = 1)$p.value
  expect_lt(abs(res$p.value - p_perm), 0.08)
  ## the one-sided default, the max, has no saddlepoint version
  expect_error(riposte_test(y ~ z | b, dat, representations = low,
                            engine = "saddlepoint", alternative = "less"), "max")
})

test_that("the asymptotic Cauchy uses one-sided normal tails", {
  d <- asymptotic_example()
  reps <- list(raw = function(v) v, rank = riposte_rank)
  res <- riposte_test(Y ~ trt | blk, d, statistic = "cauchy",
                      engine = "asymptotic", representations = reps,
                      alternative = "greater")
  sm <- riposte_score_matrix(d$Y, d$blk, reps)
  mom <- riposte_sw_moments(sm$scores, d$trt, d$blk)
  zstat <- (colSums(sm$scores * d$trt) - mom$mu) / sqrt(diag(mom$Sigma))
  p1 <- pnorm(zstat, lower.tail = FALSE)
  expect_equal(unname(res$component_p), unname(p1))
  ## the default combines the one-sided p-values with the truncated conversion;
  ## cauchy_truncation = 1 gives Liu and Xie's untruncated one
  expect_equal(res$p.value, riposte_truncated_cauchy(p1))
  res1 <- riposte_test(Y ~ trt | blk, d, statistic = "cauchy",
                       engine = "asymptotic", representations = reps,
                       alternative = "greater", cauchy_truncation = 1)
  expect_equal(res1$p.value, pcauchy(mean(tan((0.5 - p1) * pi)), lower.tail = FALSE))
})

test_that("the direction matters: harm gives a small 'less' p-value and a large 'greater' one", {
  d <- one_sided_design(3, B = 8L, nb = 20L)
  hurt <- which(d$z == 1)[seq(1, 80, by = 8)]          # 10 of 80 treated units
  y <- d$y0; y[hurt] <- y[hurt] - 4
  dat <- data.frame(y = y, z = d$z, b = d$block)
  low <- riposte_poly_reps(c(5, 14), tail = "lower")
  p_less <- riposte_test(y ~ z | b, dat, representations = low,
                         alternative = "less", nresample = 499, seed = 1)$p.value
  p_greater <- riposte_test(y ~ z | b, dat, representations = low,
                            alternative = "greater", nresample = 499, seed = 1)$p.value
  expect_lt(p_less, 0.05)
  expect_gt(p_greater, 0.5)
})

test_that("riposte_components reports one-sided mid-p values", {
  d <- one_sided_design(4)
  dat <- data.frame(y = d$y, z = d$z, b = d$block)
  reps <- list(raw = function(v) v)
  g <- riposte_components(y ~ z | b, dat, representations = reps,
                          alternative = "greater", nresample = 199, seed = 1)
  l <- riposte_components(y ~ z | b, dat, representations = reps,
                          alternative = "less", nresample = 199, seed = 1)
  ## same draws, no ties among the sums: the two one-sided mid-p values add to 1
  expect_equal(unname(g$component_p + l$component_p), 1)
  expect_output(print(g), "greater")
})

test_that("covariance adjustment accepts a one-sided alternative", {
  d <- one_sided_design(5)
  dat <- data.frame(y = d$y, z = d$z, b = d$block, x = rnorm(length(d$y)))
  res <- riposte_test(y ~ z | b, dat, adjust = ~ x, alternative = "less",
                      representations = riposte_poly_reps(c(2, 5), tail = "lower"),
                      nresample = 49, seed = 1)
  expect_equal(res$combination, "max")
  expect_true(res$p.value > 0 && res$p.value <= 1)
})

## ---- size and power (Monte Carlo; skipped on CRAN) ------------------------------

test_that("under the sharp null the one-sided max and Cauchy hold the 0.05 level", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261002)
  nsims <- 300L; alpha <- 0.05
  low <- riposte_poly_reps(c(2, 5, 14), tail = "lower")
  rej <- matrix(FALSE, nsims, 2, dimnames = list(NULL, c("max", "cauchy")))
  for (i in seq_len(nsims)) {
    d <- one_sided_design(i)
    sm <- riposte_score_matrix(d$y, d$block, low)
    mom <- riposte_sw_moments(sm$scores, d$z, d$block)
    G <- riposte_block_draws(d$z, d$block, 199L)
    rej[i, "max"] <- riposte:::riposte_max(sm$scores, d$z, d$block, moments = mom,
                                           draws = G, alternative = "less")$p.value <= alpha
    rej[i, "cauchy"] <- riposte:::riposte_cauchy(sm$scores, d$z, d$block, draws = G,
                                                 alternative = "less")$p.value <= alpha
  }
  ## SE of a rate of 0.05 over 300 experiments is about 0.0126; allow 2.5 SE
  expect_true(all(colMeans(rej) <= alpha + 2.5 * sqrt(alpha * (1 - alpha) / nsims)))
})

test_that("when some are helped and none harmed, the 'less' test rejects at most about 5 percent", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261003)
  nsims <- 300L; alpha <- 0.05
  low <- riposte_poly_reps(c(5, 14), tail = "lower")
  rej <- logical(nsims)
  for (i in seq_len(nsims)) {
    d <- one_sided_design(10000 + i)
    helped <- sample(which(d$z == 1), 6L)              # 6 of 30 treated helped by 3
    y <- d$y0; y[helped] <- y[helped] + 3
    dat <- data.frame(y = y, z = d$z, b = d$block)
    rej[i] <- riposte_test(y ~ z | b, dat, representations = low,
                           alternative = "less", nresample = 199)$p.value <= alpha
  }
  expect_lte(mean(rej), alpha + 2.5 * sqrt(alpha * (1 - alpha) / nsims))
})

test_that("against harm to a few, the one-sided bottom-weighted test beats the two-sided difference in means", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261004)
  nsims <- 100L; alpha <- 0.05
  low <- riposte_poly_reps(c(5, 14, 41), tail = "lower")
  raw <- list(raw = function(v) v)
  rej <- matrix(FALSE, nsims, 2, dimnames = list(NULL, c("bottom", "means")))
  for (i in seq_len(nsims)) {
    d <- one_sided_design(20000 + i, B = 20L, nb = 40L)
    hurt <- sample(which(d$z == 1), 20L)               # 20 of 400 treated, by 3 SD
    y <- d$y0; y[hurt] <- y[hurt] - 3
    dat <- data.frame(y = y, z = d$z, b = d$block)
    rej[i, "bottom"] <- riposte_test(y ~ z | b, dat, representations = low,
                                     alternative = "less",
                                     nresample = 199)$p.value <= alpha
    rej[i, "means"] <- riposte_test(y ~ z | b, dat, representations = raw,
                                    statistic = "quadratic",
                                    engine = "asymptotic")$p.value <= alpha
  }
  rr <- colMeans(rej)
  ## In the vignette's version of this design (100 pilots) the one-sided
  ## bottom-weighted test found the harm in 91 and a one-sided test of the
  ## average in 57; the two-sided test of the average does worse still.
  expect_gt(rr[["bottom"]], rr[["means"]] + 0.15)
})

## ---- agreement with CMRSS ---------------------------------------------------------

test_that("the one-sided max of bottom-weighted scores is CMRSS's test that no treated unit was harmed", {
  ## CMRSS (Kim, Su, Bowers, and Li) tests hypotheses about how many treated
  ## units had effects above a threshold. Its test that the largest effect among
  ## treated units is at most 0 (k = number treated, c = 0), applied to -y with
  ## top-weighted polynomial scores, is the hypothesis that no treated unit was
  ## harmed. With blocks of equal size its statistic is the largest of the
  ## standardized treated-score sums, the same statistic riposte's one-sided max
  ## computes from bottom-weighted scores of y with alternative = "less".
  skip_if_not_installed("CMRSS")
  skip_if_not_installed("highs")
  d <- one_sided_design(8, B = 6L, nb = 20L)
  hurt <- which(d$z == 1)[c(2, 17, 40)]
  y <- d$y0; y[hurt] <- y[hurt] - 3
  zetas <- c(5, 14)
  rip <- riposte_test(y ~ z | b, data.frame(y = y, z = d$z, b = d$block),
                      representations = riposte_poly_reps(zetas, tail = "lower"),
                      alternative = "less", nresample = 99, seed = 1)
  methods <- lapply(zetas, function(r) lapply(seq_len(6), function(i)
    list(name = "Polynomial", r = r, std = TRUE, scale = FALSE)))
  cm <- suppressWarnings(CMRSS::pval_comb_block(d$z, -y, k = sum(d$z), c = 0,
                                                block = d$block,
                                                methods.list.all = methods,
                                                null.max = 200,
                                                opt.method = "ILP_highs"))
  expect_equal(rip$statistic, unname(cm[["test.stat"]]), tolerance = 1e-8)
})
