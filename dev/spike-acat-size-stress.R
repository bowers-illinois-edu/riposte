## spike-acat-size-stress.R --- stress the unadjusted analytic ACAT to find a
## false-positive regime the 40-block check might have missed, and pin the driver.
##
## The first diagnostic (spike-acat-size.R) found only mild anti-conservatism in
## the unadjusted case at 40 blocks, and showed dependence alone does not break
## the analytic tail. Two follow-ups here:
##   PART A  few-blocks sweep on the FULL riposte test, averaged over outcome
##           draws -- the finite-sample regime where the Cauchy approximation is
##           furthest from its asymptotics and the permutation marginals are
##           coarsest. Analytic vs permutation (the size-valid reference).
##   PART B  vary REDUNDANCY directly (cheap Gaussian copula, exactly-uniform
##           marginals): k near-duplicate components out of d. Tests whether the
##           collinear representations (raw vs tanh ~ 0.98) drive the small
##           Rung-2b inflation. Theory predicts NOT: the average of perfectly
##           dependent standard Cauchys is again standard Cauchy, so strong
##           redundancy should preserve calibration. Part B checks that prediction.

source("dev/acat-helpers.R")
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625L)
ALPHA <- c(0.10, 0.05, 0.01)

## ---- PART A: few-blocks sweep, averaged over outcome draws -------------------
## Fixed bs/m; vary the number of blocks. For each block count, draw several
## outcomes y (so the result is not conditional on one realized correlation
## structure) and, per y, several fresh observed assignments z under H0.
run_avg <- function(n_blocks, bs = 6L, m = 3L, y_reps = 50L, z_reps = 100L, B = 999L) {
  N <- y_reps * z_reps
  p_an <- numeric(N); p_pm <- numeric(N); k <- 0L
  for (yr in seq_len(y_reps)) {
    dn <- mk_design(n_blocks, bs, m)                       # fresh outcome
    scores <- riposte_score_matrix(dn$y, dn$block)$scores
    idx <- split(seq_len(dn$n), dn$block)
    for (zr in seq_len(z_reps)) {
      z <- numeric(dn$n); for (ix in idx) z[ix[sample.int(length(ix), m)]] <- 1
      res <- riposte_cauchy(scores, z, dn$block, nresample = B)
      k <- k + 1L
      p_pm[k] <- res$p.value
      p_an[k] <- analytic_acat_p(matrix(res$component_midp, 1))
    }
  }
  list(p_an = p_an, p_pm = p_pm, N = N)
}

cat("=== PART A: few-blocks stress (full riposte, averaged over outcomes) ===\n")
cat("  bs=6, m=3; analytic vs PERMUTATION (size-valid) FPR; .10/.05/.01\n")
cat(sprintf("  %-10s %6s %26s %26s\n", "n_blocks", "N", "analytic", "permutation"))
for (nb in c(2L, 4L, 8L, 16L, 40L)) {
  r <- run_avg(nb)
  cat(sprintf("  %-10d %6d   %.3f/%.3f/%.3f          %.3f/%.3f/%.3f\n",
              nb, r$N,
              fpr(r$p_an, .10), fpr(r$p_an, .05), fpr(r$p_an, .01),
              fpr(r$p_pm, .10), fpr(r$p_pm, .05), fpr(r$p_pm, .01)))
}
cat(sprintf("  (MC SE at .05 ~ %.4f for N=5000)\n", fpr_se(.05, 5000)))

## ---- PART B: redundancy sweep (cheap, exactly-uniform marginals) -------------
## d components; the first k share a near-common factor (pairwise corr ~0.99),
## the rest independent. k = 1 -> all independent; k = d -> all near-duplicate.
cat("\n=== PART B: does REDUNDANCY inflate the analytic ACAT? ===\n")
cat("  d=8 components, k near-duplicates (corr ~0.99), exactly-uniform marginals\n")
Nbig <- 50000L; d <- 8L; rho_dup <- 0.995
cat(sprintf("  %-8s %26s\n", "k dup", "analytic FPR .10/.05/.01"))
for (k in c(1L, 2L, 4L, 6L, 8L)) {
  f1 <- rnorm(Nbig)
  Z <- matrix(rnorm(Nbig * d), Nbig, d)                    # independent base
  if (k >= 1L) for (j in seq_len(k))                       # tie first k to f1
    Z[, j] <- rho_dup * f1 + sqrt(1 - rho_dup^2) * Z[, j]
  p_an <- analytic_acat_p(unif_marginals_from_Z(Z))
  cat(sprintf("  %-8d %.4f / %.4f / %.4f\n", k,
              fpr(p_an, .10), fpr(p_an, .05), fpr(p_an, .01)))
}
cat(sprintf("  (MC SE at .05 ~ %.4f for N=5e4)\n", fpr_se(.05, Nbig)))
