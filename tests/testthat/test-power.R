## test-power.R
##
## POWER DIRECTIONS (spec section 8). The point of combining representations is to
## detect effects a difference of means misses, and the point of the quadratic
## (which uses the full permutation covariance) over the Cauchy (which uses only
## the marginals) is to pool evidence that is spread across correlated
## representations. These tests confirm the qualitative orderings the proposal's
## Table reports, on a smaller design so they run quickly:
##
##   - Cancelling effect (half the blocks +delta, half -delta; the average
##     cancels): the QUADRATIC detects it, the CAUCHY does NOT, and a difference
##     of means is blind. The signal is individually weak in each distance
##     representation but correlated across them, so only the covariance-using
##     quadratic recovers it. This is the package's signature claim.
##   - Scale effect (treatment changes spread, not mean): the quadratic detects
##     it; a difference of means is blind.
##   - Sparse upper-tail effect (a few units helped a lot): a combined test beats
##     a single rank (Wilcoxon) test.
##   - Covariates of pure noise: covariance adjustment, refit inside every
##     re-randomization, rejects about as often as the unadjusted test. The main
##     vignette quotes the two rates this test produces at its seed, 0.85 and 0.83.
##
## These are Monte-Carlo rejection rates, so they are slow and skipped on CRAN.
## How the thresholds are set: each is the calibrated rejection rate minus a
## generous margin (several Monte-Carlo SEs), so the test holds across seeds and
## platforms but still fails if a combination loses the property it is supposed to
## have. The Monte-Carlo SE of a rate p over n sims is sqrt(p (1 - p) / n). The
## calibrated rate and SE are given next to each assertion so the margin is
## visible and reproducible (measured 2026-06-23 with the seeds below).

## all four p-values from one shared set of draws: the two combinations, plus the
## raw and rank marginals (a stratified difference of means and a single rank
## test, both referred to the same randomization distribution)
power_pvalues <- function(y, z, block, nres) {
  sm <- riposte_score_matrix(y, block)
  mom <- riposte_sw_moments(sm$scores, z, block)
  G <- riposte_block_draws(z, block, nres)
  Tmat <- crossprod(G, sm$scores)
  marg <- apply(Tmat, 2, function(col) mean(abs(col) >= abs(col[1])))
  c(quadratic = riposte_quadratic(sm$scores, z, block, moments = mom, draws = G)$p.value,
    cauchy    = riposte_cauchy(sm$scores, z, block, draws = G)$p.value,
    diffmeans = marg[["raw"]],
    wilcox    = marg[["rank"]])
}

assign_within_block <- function(block) {
  z <- integer(length(block))
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(length(ix), length(ix) %/% 2L)] <- 1L
  }
  z
}

reject_rates <- function(make_data, nsims, nres, alpha = 0.05) {
  P <- t(replicate(nsims, {
    d <- make_data()
    power_pvalues(d$y, d$z, d$block, nres)
  }))
  colMeans(P <= alpha)
}

test_that("the quadratic detects a cancelling effect that the Cauchy and diff of means miss", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260623)
  B <- 24L; nb <- 50L; delta <- 1.3
  block <- factor(rep(seq_len(B), each = nb))
  make <- function() {
    z <- assign_within_block(block)
    y0 <- rnorm(B * nb)
    pos <- as.integer(block) <= (B / 2)               # first half +delta, rest -delta
    list(y = y0 + z * ifelse(pos, delta, -delta), z = z, block = block)
  }
  rr <- reject_rates(make, nsims = 100L, nres = 199L)

  ## the quadratic recovers the cancelling effect; the Cauchy and the difference
  ## of means do not. Calibrated rates (n = 100, SE ~ sqrt(p(1-p)/100)):
  ##   quadratic 0.44 (SE 0.050), cauchy 0.04 (SE 0.020), diffmeans 0.00.
  expect_gt(rr[["quadratic"]], 0.25)                  # 0.44 - ~4 SE
  expect_lt(rr[["cauchy"]], 0.20)                     # 0.04 + ~8 SE
  expect_lt(rr[["diffmeans"]], 0.10)                  # average cancels: blind
  ## the covariance-using quadratic beats the marginal-only Cauchy by a clear
  ## margin (calibrated gap ~0.40)
  expect_gt(rr[["quadratic"]] - rr[["cauchy"]], 0.15)
})

test_that("the quadratic detects a scale change the difference of means misses", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260624)
  B <- 16L; nb <- 25L; sigma1 <- 1.5
  block <- factor(rep(seq_len(B), each = nb))
  make <- function() {
    z <- assign_within_block(block)
    y0 <- rnorm(B * nb); y1 <- rnorm(B * nb, 0, sigma1)  # spread changes, mean does not
    list(y = ifelse(z == 1, y1, y0), z = z, block = block)
  }
  rr <- reject_rates(make, nsims = 80L, nres = 199L)
  ## Calibrated rates (n = 80): quadratic 0.99 (SE 0.011), diffmeans 0.05 (SE 0.024).
  expect_gt(rr[["quadratic"]], 0.85)                  # 0.99 - ~13 SE (loose on purpose)
  expect_lt(rr[["diffmeans"]], 0.15)                  # diff of means blind to a scale change
})

test_that("a combined test beats a single rank test on a sparse upper-tail effect", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260625)
  B <- 24L; nb <- 50L; frac <- 0.05; tau <- 2.5
  block <- factor(rep(seq_len(B), each = nb))
  make <- function() {
    z <- assign_within_block(block)
    y0 <- rnorm(B * nb)
    ## a few units per block (ceiling(frac * nb) = 3 of 50) helped a lot
    resp <- ave(seq_along(block), block,
                FUN = function(ii) as.integer(seq_along(ii) <= ceiling(frac * length(ii))))
    list(y = y0 + z * resp * tau, z = z, block = block)
  }
  rr <- reject_rates(make, nsims = 80L, nres = 199L)
  ## Calibrated rates (n = 80): quadratic 0.85 (SE 0.040), wilcox 0.34 (SE 0.053).
  expect_gt(rr[["quadratic"]], 0.50)                  # 0.85 - ~9 SE
  ## the combined quadratic beats the single rank (Wilcoxon) test (calibrated
  ## gap ~0.51)
  expect_gt(rr[["quadratic"]] - rr[["wilcox"]], 0.15)
})

test_that("adjusting for covariates of pure noise loses little power", {
  skip_on_cran()
  ## The main vignette's design (12 blocks of 30, half treated, treatment widens
  ## the spread) with a smaller effect, sd 1 -> 1.35, so that power is below 1
  ## and a loss would show. x1 and x2 are unrelated to the outcome, so the
  ## controls-only ridge fit has nothing to remove; any power it loses comes from
  ## fitting noise. Both tests see the same data in every replication, so the
  ## difference in rejection rates is a paired comparison.
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261005)
  B <- 12L; nb <- 30L; nsims <- 100L
  block <- factor(rep(seq_len(B), each = nb))
  rej <- t(replicate(nsims, {
    z <- assign_within_block(block)
    y <- ifelse(z == 1, rnorm(B * nb, sd = 1.35), rnorm(B * nb))
    d <- data.frame(Y = y, trt = z, blk = block,
                    x1 = rnorm(B * nb), x2 = rnorm(B * nb))
    c(unadj = riposte_test(Y ~ trt | blk, d, nresample = 99)$p.value <= 0.05,
      adj   = riposte_test(Y ~ trt | blk, d, adjust = ~ x1 + x2,
                           nresample = 99)$p.value <= 0.05)
  }))
  rr <- colMeans(rej)
  ## Calibrated rates (2026-10-02, n = 400, set.seed(1)): unadjusted 0.845,
  ## adjusted 0.835. At this seed (n = 100): 0.85 and 0.83. Over seeds 2, 3, 4
  ## (n = 100 each) the unadjusted rate exceeded the adjusted by 0.03, 0.05, 0.00.
  ## The tests disagree in about 5 to 8 percent of replications, so the SE of the
  ## paired difference at n = 100 is about sqrt(0.07 / 100) = 0.026.
  expect_gt(rr[["adj"]], 0.65)                         # 0.83 - ~5 SE (SE 0.038)
  expect_lt(rr[["unadj"]] - rr[["adj"]], 0.10)         # 0.02 + ~3 SE of the difference
})

