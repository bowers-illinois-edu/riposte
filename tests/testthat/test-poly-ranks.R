## test-poly-ranks.R
##
## The statistical content:
##   - the polynomial rank score of a unit at within-block rank k in a block of
##     n_b units is (k / (n_b + 1))^(zeta - 1), ties at their average rank, as in
##     Kim, Li, and Bowers and in Bowers and Burton's rep_poly_rank();
##   - unlike Stephenson's choose(k - 1, zeta - 1), the score stays in [0, 1) in
##     every block, so blocks of different sizes enter on the same scale --- the
##     reason the paper prefers this form;
##   - zeta = 2 is the rank rescaled: after within-block centring it is
##     proportional to the rank representation when blocks share a size, so it
##     duplicates `rank` in riposte_reps_default();
##   - for a block large enough that rank / (n_b + 1) is close to uniform on
##     (0, 1), the correlation of the scores at zeta and zeta' is close to
##     2 sqrt(u v) / (u + v) with u = 2 zeta - 1 and v = 2 zeta' - 1, which
##     depends only on the ratio v / u. That is the fact the vignette uses to
##     space a user's zeta values;
##   - the set flows through riposte_test with any zeta the user supplies, and
##     the combined test holds size under the sharp null;
##   - on a spread change with equal means, the top-weighted scores detect what
##     the rank (zeta = 2) misses.

test_that("the polynomial score equals (rank / (n + 1))^(zeta - 1), ties averaged", {
  y <- c(3.1, -0.4, 2.2, 2.2, 7.5, 0.0)
  rk <- rank(y, ties.method = "average")              # the two 2.2s share 3.5
  n <- length(y)
  expect_equal(riposte_poly_scores(y, 2), rk / (n + 1))
  expect_equal(riposte_poly_scores(y, 7), (rk / (n + 1))^6)
  ## zeta need not be a whole number: the power is defined for any zeta >= 1,
  ## unlike choose() in the Stephenson score
  expect_equal(riposte_poly_scores(y, 2.5), (rk / (n + 1))^1.5)
  ## zeta = 1 gives the constant 1, which carries no information
  expect_equal(riposte_poly_scores(y, 1), rep(1, n))
  expect_error(riposte_poly_scores(y, 0.5), "zeta")
  expect_error(riposte_poly_scores(y, c(2, 3)), "zeta")
  expect_error(riposte_poly_scores(y, NA), "zeta")
})

test_that("the polynomial score stays in [0, 1) whatever the block size", {
  ## the Stephenson score at zeta = 22 reaches choose(49, 21), about 3.9e13, in a
  ## block of 50; the polynomial score never reaches 1 because rank <= n < n + 1
  for (n in c(5L, 50L, 500L)) {
    y <- rnorm(n)
    for (zeta in c(2, 7, 22, 60)) {
      s <- riposte_poly_scores(y, zeta)
      expect_true(all(s >= 0 & s < 1))
    }
  }
  expect_gt(choose(49, 21), 1e13)
})

test_that("zeta = 2 duplicates the rank after centring when blocks share a size", {
  set.seed(11)
  y <- rnorm(40)
  equal <- factor(rep(1:4, each = 10))
  reps <- list(rank = riposte_rank, poly2 = function(y) riposte_poly_scores(y, 2))
  cs <- riposte_centred_scores(y, equal, reps)
  ## every block has n_b = 10, so the centred poly2 column is the centred rank
  ## column times 1 / 11 in every block
  expect_equal(cs[, "poly2"], cs[, "rank"] / 11)

  ## with unequal block sizes the factor 1 / (n_b + 1) differs by block, so the
  ## two columns are no longer exactly proportional
  unequal <- factor(rep(1:2, times = c(10, 30)))
  cu <- riposte_centred_scores(y, unequal, reps)
  ratio <- cu[, "poly2"] / cu[, "rank"]
  ratio <- ratio[is.finite(ratio)]
  expect_gt(diff(range(ratio)), 0.01)
})

test_that("the correlation between two zetas follows 2 sqrt(uv) / (u + v)", {
  ## For one block, the permutation correlation of the two treated-unit sums
  ## equals the correlation of the two centred score vectors, and with no ties
  ## that depends only on n_b. Treating rank / (n_b + 1) as uniform on (0, 1),
  ## corr(X^a, X^b) = sqrt((2a + 1)(2b + 1)) / (a + b + 1) with a = zeta - 1.
  closed <- function(z1, z2) {
    u <- 2 * z1 - 1; v <- 2 * z2 - 1
    2 * sqrt(u * v) / (u + v)
  }
  y <- seq_len(30)                                     # ranks 1..30, no ties
  pairs <- list(c(2, 7), c(7, 12), c(17, 22), c(2, 5), c(5, 14))
  for (p in pairs) {
    exact <- cor(riposte_poly_scores(y, p[1]), riposte_poly_scores(y, p[2]))
    expect_equal(exact, closed(p[1], p[2]), tolerance = 0.01)
  }
  ## equal ratios of 2 zeta - 1 give equal neighbour correlations: 3, 9, 27 are
  ## 2 zeta - 1 at zeta = 2, 5, 14
  expect_equal(closed(2, 5), closed(5, 14))
  ## equal steps in zeta do not: the correlation of neighbours climbs toward 1
  expect_lt(closed(2, 7), closed(7, 12))
  expect_lt(closed(7, 12), closed(17, 22))
})

test_that("riposte_poly_reps() builds a named set from any zeta the user supplies", {
  reps <- riposte_poly_reps(c(3, 9, 20))
  expect_identical(names(reps), c("poly3", "poly9", "poly20"))
  expect_true(all(vapply(reps, is.function, logical(1))))
  y <- c(4, 1, 8, 2, 9)
  expect_equal(reps$poly9(y), riposte_poly_scores(y, 9))

  set.seed(3)
  block <- factor(rep(1:5, each = 8)); N <- 40
  z <- integer(N)
  for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(8, 4)] <- 1L }
  d <- data.frame(Y = rnorm(N) + 2 * z * (seq_len(N) %% 8 == 1), trt = z, blk = block)
  r <- riposte_test(Y ~ trt | blk, d, representations = reps,
                    statistic = "quadratic", nresample = 199, seed = 1)
  expect_s3_class(r, "riposte_test")
  expect_setequal(r$kept, names(reps))

  ## zeta = 1 is constant within every block, so riposte drops it
  r1 <- riposte_test(Y ~ trt | blk, d, representations = riposte_poly_reps(c(1, 6)),
                     statistic = "quadratic", nresample = 99, seed = 1)
  expect_identical(r1$dropped, "poly1")

  ## the set adds to the six defaults with c()
  both <- c(riposte_reps_default(), riposte_poly_reps(c(7, 12)))
  expect_length(both, 8L)
})

test_that("the combined polynomial-rank test holds size under the sharp null", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261002)
  B <- 6L; nb <- 12L; nsims <- 150L
  block <- factor(rep(seq_len(B), each = nb))
  reps <- riposte_poly_reps(c(2, 5, 11))
  rej <- numeric(nsims)
  for (i in seq_len(nsims)) {
    y <- rnorm(B * nb)                                  # sharp null
    z <- integer(B * nb)
    for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L }
    d <- data.frame(Y = y, trt = z, blk = block)
    rej[i] <- riposte_test(Y ~ trt | blk, d, representations = reps,
                           statistic = "quadratic", nresample = 149)$p.value <= 0.05
  }
  rate <- mean(rej)
  expect_lt(rate, 0.10)
  expect_gt(rate, 0.01)
})

test_that("top-weighted scores detect a spread change the rank misses", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261003)
  ## the vignette's design: 12 blocks of 30, half treated, treatment raises the
  ## standard deviation from 1 to 1.8 and leaves the mean at 0.
  ## Calibrated (n = 200 at this seed): zeta = 2 alone 0.045, zeta = 5 and 14 0.99.
  B <- 12L; nb <- 30L; nsims <- 60L
  block <- factor(rep(seq_len(B), each = nb))
  power <- function(reps) {
    p <- numeric(nsims)
    for (i in seq_len(nsims)) {
      z <- integer(B * nb)
      for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L }
      y <- ifelse(z == 1, rnorm(B * nb, sd = 1.8), rnorm(B * nb))
      d <- data.frame(Y = y, trt = z, blk = block)
      p[i] <- riposte_test(Y ~ trt | blk, d, representations = reps,
                           statistic = "quadratic", nresample = 199)$p.value <= 0.05
    }
    mean(p)
  }
  rank_only <- power(riposte_poly_reps(2))
  top <- power(riposte_poly_reps(c(5, 14)))
  expect_lt(rank_only, 0.15)                           # a spread change leaves ranks' mean alone
  expect_gt(top, 0.85)
})
