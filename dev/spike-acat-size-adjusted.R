## spike-acat-size-adjusted.R --- the covariance-adjustment FPR diagnostic.
##
## THE QUESTION. The remembered false-positive-rate problem is most likely the
## ADJUSTED case (refit-per-permutation), not the unadjusted one. This script,
## under the sharp null with controls-only covariance adjustment, computes four
## combined p-values per dataset and asks:
##   (1) does the documented metric failure reproduce -- draws-only mu/Sigma making
##       the QUADRATIC anti-conservative -- and does the POOLED metric restore it?
##       (the "is the symmetric-pooled-metric the fix" part)
##   (2) does the analytic Liu-Xie ACAT tail break size under adjustment, vs the
##       permutation calibration of the same residual ACAT statistic? (the "does
##       analytic ACAT break there" part)
## The Cauchy combination uses marginal mid-p, NOT mu/Sigma, so a metric failure
## cannot hit it; the marginals stay uniform under refit + pooled exchangeability,
## so the prediction is that analytic ACAT behaves much as it did unadjusted. We
## test that rather than assume it.
##
## Reproduced number to watch: the theory note reports draws-only metric -> size
## ~0.12 at SMALL nresample; pooled -> ~0.05. We sweep B to see that O(1/B) effect.

source("dev/acat-helpers.R")
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625L)
ALPHA <- c(0.10, 0.05, 0.01)

## one adjusted, block-randomized null dataset: covariates predict y, NO effect
n_blocks <- 20L; bs <- 6L; m <- 3L; d_cov <- 8L
n <- n_blocks * bs
block <- factor(rep(seq_len(n_blocks), each = bs))
idx <- split(seq_len(n), block)
lam0 <- 1                                                   # fixed ridge penalty (size is lambda-invariant)
learner <- riposte_ridge_learner(lambda = lam0)
draw_z <- function() { z <- numeric(n); for (ix in idx) z[ix[sample.int(length(ix), m)]] <- 1; z }
mk_data <- function() {
  X <- matrix(rnorm(n * d_cov), n, d_cov)
  y <- as.numeric(X %*% rep(0.7, d_cov)) + rnorm(n)        # predictive covariates, sharp null
  list(X = X, y = y)
}

## one dataset -> four combined p-values from the refit-per-permutation residuals
four_pvalues <- function(B) {
  dat <- mk_data(); z <- draw_z()
  draws <- riposte_block_draws(z, block, nresample = B)    # col 1 = observed
  Tm <- riposte_residual_perm_stats(dat$y, dat$X, block, draws, learner)$stats
  ## quadratic with the POOLED metric (riposte's exact default) ...
  q_pool <- riposte_quadratic_from_T(Tm, colMeans(Tm), stats::cov(Tm))$p.value
  ## ... and with the DRAWS-ONLY metric (the documented anti-conservative bug)
  Td <- Tm[-1, , drop = FALSE]
  q_draw <- riposte_quadratic_from_T(Tm, colMeans(Td), stats::cov(Td))$p.value
  ## Cauchy: permutation calibration (riposte default) vs analytic Liu-Xie tail
  cau <- riposte_cauchy_from_T(Tm)
  c(q_pool = q_pool, q_draw = q_draw,
    cau_perm = cau$p.value, cau_anal = analytic_acat_p(matrix(cau$component_midp, 1)))
}

cat("=== Covariance-adjustment FPR diagnostic (sharp null, refit-per-permutation) ===\n")
cat(sprintf("  %d blocks x %d (m=%d), %d covariates, ridge lambda=%g\n",
            n_blocks, bs, m, d_cov, lam0))
run_B <- function(B, Nsim) {
  P <- matrix(NA_real_, Nsim, 4, dimnames = list(NULL, c("q_pool","q_draw","cau_perm","cau_anal")))
  for (i in seq_len(Nsim)) P[i, ] <- four_pvalues(B)
  P
}
for (cfg in list(c(B = 49L, Nsim = 2000L), c(B = 199L, Nsim = 2000L))) {
  B <- cfg["B"]; Nsim <- cfg["Nsim"]
  P <- run_B(B, Nsim)
  cat(sprintf("\n[B=%d, Nsim=%d; MC SE at .05 ~ %.4f]\n", B, Nsim, fpr_se(.05, Nsim)))
  cat(sprintf("  %-26s  %s\n", "combination", "FPR .10/.05/.01"))
  labs <- c(q_pool = "quadratic POOLED (exact)", q_draw = "quadratic draws-only (bug)",
            cau_perm = "cauchy permutation (exact)", cau_anal = "cauchy ANALYTIC (Liu-Xie)")
  for (k in colnames(P))
    cat(sprintf("  %-26s  %.3f / %.3f / %.3f\n", labs[k],
                fpr(P[, k], .10), fpr(P[, k], .05), fpr(P[, k], .01)))
}
