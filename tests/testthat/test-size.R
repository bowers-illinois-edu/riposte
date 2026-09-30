## test-size.R
##
## THE CENTRAL CLAIM. Under the sharp null of no effect, every combination must
## hold the nominal 0.05 level. This is the property the whole package rests on:
## the p-values come from the randomization distribution, so they are valid in
## finite samples by construction. We check it by simulation --- draw data with no
## effect, randomize, test, and confirm each combination rejects at about 0.05.
##
## This is a Monte-Carlo check, so it is slow and skipped on CRAN. The rejection
## rate is itself estimated with error ~ sqrt(0.05 * 0.95 / nsims); the tolerance
## band below is wide enough to pass under the null and narrow enough to catch a
## combination that does not hold its size (e.g. an anti-conservative quadratic
## from Monte-Carlo moments, the bug the closed-form moments fix).

test_that("under the sharp null every combination holds the 0.05 level", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260623)
  B <- 6L; nb <- 10L; nsims <- 300L; nres <- 199L; alpha <- 0.05
  block <- factor(rep(seq_len(B), each = nb))

  rej <- matrix(0, nsims, 3,
                dimnames = list(NULL, c("quadratic", "cauchy", "max")))
  for (i in seq_len(nsims)) {
    y <- rnorm(B * nb)                 # sharp null: no effect anywhere
    z <- integer(B * nb)
    for (b in levels(block)) {
      ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L
    }
    sm <- riposte_score_matrix(y, block)
    mom <- riposte_sw_moments(sm$scores, z, block)
    G <- riposte_block_draws(z, block, nres)     # one shared set of draws
    rej[i, "quadratic"] <- riposte_quadratic(sm$scores, z, block,
                                             moments = mom, draws = G)$p.value <= alpha
    rej[i, "cauchy"] <- riposte_cauchy(sm$scores, z, block, draws = G)$p.value <= alpha
    rej[i, "max"] <- riposte_max(sm$scores, z, block,
                                 moments = mom, draws = G)$p.value <= alpha
  }
  rr <- colMeans(rej)
  ## holds 0.05: rejection rate near nominal, not inflated. Calibrated rates
  ## (n = 300, SE ~ sqrt(0.05 * 0.95 / 300) = 0.0126): quadratic 0.050, cauchy
  ## 0.053, max 0.060. The band 0.015-0.09 is about +/- 3 SE around 0.05 -- wide
  ## enough to pass under the null, narrow enough to catch an anti-conservative
  ## combination (e.g. the Monte-Carlo-moments quadratic the closed form fixes).
  for (nm in colnames(rej)) {
    expect_lt(rr[[nm]], 0.09)
    expect_gt(rr[[nm]], 0.015)
  }
})

test_that("the screen (the package default) holds the level end to end", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260626)
  B <- 6L; nb <- 10L; nsims <- 200L
  block <- factor(rep(seq_len(B), each = nb))
  rej <- numeric(nsims)
  for (i in seq_len(nsims)) {
    y <- rnorm(B * nb)
    z <- integer(B * nb)
    for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L }
    d <- data.frame(Y = y, trt = z, blk = block)
    ## statistic = "screen" is the default: this exercises the shrink-then-quadratic
    ## path that the per-combination tests above do not
    rej[i] <- riposte_test(Y ~ trt | blk, d, statistic = "screen", nresample = 149)$p.value <= 0.05
  }
  rate <- mean(rej)
  expect_lt(rate, 0.10)
  expect_gt(rate, 0.01)
})
