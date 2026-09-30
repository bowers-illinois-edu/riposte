## test-stephenson.R
##
## The statistical content:
##   - the Stephenson score is choose(rank - 1, r - 1); r = 2 is the Wilcoxon
##     (rank - 1) score, larger r upweights the top ranks;
##   - the scores are an open representation set that flows through riposte_test;
##   - the combined Stephenson test holds size under the sharp null;
##   - on a SPARSE UPPER-TAIL effect (a few units helped a lot) the combined
##     Stephenson test (several r) beats a single Wilcoxon (r = 2) -- the reason
##     to combine tail weightings instead of committing to one.

test_that("Stephenson scores equal choose(rank - 1, r - 1)", {
  set.seed(1)
  y <- rnorm(10)
  rk <- rank(y, ties.method = "average")
  expect_equal(riposte_stephenson_scores(y, 2), choose(rk - 1, 1))   # Wilcoxon = rank - 1
  expect_equal(riposte_stephenson_scores(y, 6), choose(rk - 1, 5))
  ## larger r concentrates on fewer top units: a unit scores 0 until its rank
  ## reaches r, so the number of nonzero scores shrinks as r grows
  nz <- function(r) sum(riposte_stephenson_scores(y, r) > 0)
  expect_true(nz(2) > nz(6) && nz(6) > nz(10))
  expect_error(riposte_stephenson_scores(y, 0), "r` must be")
  ## a fractional r would silently round inside choose(); reject it instead
  expect_error(riposte_stephenson_scores(y, 2.5), "whole number")
})

test_that("Stephenson reps are a named representation set that riposte_test accepts", {
  reps <- riposte_stephenson_reps(c(2, 6, 10))
  expect_setequal(names(reps), c("S2", "S6", "S10"))
  expect_true(all(vapply(reps, is.function, logical(1))))

  set.seed(2)
  block <- factor(rep(1:5, each = 8)); N <- 40
  z <- integer(N); for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(8, 4)] <- 1L }
  d <- data.frame(Y = rnorm(N) + 2 * z * (seq_len(N) %% 8 == 1), trt = z, blk = block)
  r <- riposte_test(Y ~ trt | blk, d, representations = reps, statistic = "quadratic",
                    nresample = 199, seed = 1)
  expect_s3_class(r, "riposte_test")
  expect_setequal(c(r$kept, r$dropped), names(reps))
  expect_true(r$p.value > 0 && r$p.value <= 1)
})

test_that("the combined Stephenson test holds size under the sharp null", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260624)
  B <- 6L; nb <- 12L; nsims <- 150L
  block <- factor(rep(seq_len(B), each = nb))
  reps <- riposte_stephenson_reps(c(2, 6, 10))
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

test_that("combined Stephenson beats a single Wilcoxon on a sparse upper-tail effect", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260625)
  ## Calibrated (n = 80): single Wilcoxon (S2) 0.34, combined Stephenson 0.81.
  B <- 24L; nb <- 50L; frac <- 0.05; tau <- 2.5
  block <- factor(rep(seq_len(B), each = nb))
  power <- function(nsims, reps, stat) {
    p <- numeric(nsims)
    for (i in seq_len(nsims)) {
      z <- integer(B * nb)
      for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L }
      resp <- ave(seq_along(block), block,
                  FUN = function(ii) as.integer(seq_along(ii) <= ceiling(frac * length(ii))))
      y <- rnorm(B * nb) + z * resp * tau
      d <- data.frame(Y = y, trt = z, blk = block)
      p[i] <- riposte_test(Y ~ trt | blk, d, representations = reps,
                           statistic = stat, nresample = 199)$p.value <= 0.05
    }
    mean(p)
  }
  wilcox <- power(70L, riposte_stephenson_reps(2), "max")
  combined <- power(70L, riposte_stephenson_reps(c(2, 6, 10)), "quadratic")
  expect_gt(combined, 0.55)                            # calibrated 0.81
  expect_gt(combined - wilcox, 0.20)                   # calibrated gap ~0.47
})
