## spike-acat-fast.R --- the fast, draws-free Cauchy combination for the
## unadjusted test: SPA per-representation marginals -> Liu-Xie analytic tail.
##
## THE POINT. riposte's Cauchy combination needs each representation's two-sided
## permutation p-value, then combines them. The marginals are the expensive part
## (they need the within-block draws). The single-representation saddlepoint
## (dev/saddlepoint-core.R) returns each marginal in closed form, with NO draws;
## the Liu-Xie analytic tail then combines them in closed form. The size diagnostic
## (dev/spike-acat-size*.R) showed this analytic combination controls the FPR, so
## the whole Cauchy p-value can be computed without re-randomizing at all. This
## script validates the fast path against riposte's brute-force permutation Cauchy
## on the 44-block / N=2,200 example, separating the two approximation sources:
##   - SPA marginal error: SPA two-sided p_j vs brute-force two-sided mid-p
##   - analytic-vs-permutation calibration: analytic ACAT vs permutation of T_c
## and times the draws-free path against the brute-force one.

source("dev/saddlepoint-core.R")
source("dev/acat-helpers.R")
suppressMessages(devtools::load_all(".", quiet = TRUE))

## 44 blocks of 50, 25 treated, modest effect so the marginals are informative
set.seed(20260625L)
n_blocks <- 44L; bs <- 50L; m <- 25L; n <- n_blocks * bs
block <- factor(rep(seq_len(n_blocks), each = bs))
z <- numeric(n); for (ix in split(seq_len(n), block)) z[ix[seq_len(m)]] <- 1
y <- rnorm(n) + 0.08 * z

S <- riposte_score_matrix(y, block)                        # default 6 representations
scores <- S$scores; reps <- colnames(scores)

## ---- fast path: SPA two-sided marginal per representation, no draws ----------
tfast <- system.time({
  p_spa <- vapply(seq_along(reps), function(j) {
    spa <- make_spa(scores[, j], z, block)
    spa$two_sided(sum(z * scores[, j]))                    # observed two-sided p
  }, numeric(1))
  p_fast <- analytic_acat_p(matrix(p_spa, 1))              # Liu-Xie analytic combine
})[["elapsed"]]
names(p_spa) <- reps

## ---- reference: brute-force permutation marginals and permutation Cauchy ------
B <- 100000L
tbrute <- system.time({
  draws <- riposte_block_draws(z, block, nresample = B)
  Tmat  <- crossprod(draws, scores)
  colnames(Tmat) <- reps
  ## riposte's two-sided mid-p marginal at the observed assignment (row 1)
  p_mid <- apply(Tmat, 2, function(col) riposte_midp(col)[1])
  ## riposte's permutation Cauchy p-value (calibrates T_c by permutation)
  cau   <- riposte_cauchy_from_T(Tmat)
})[["elapsed"]]
## analytic ACAT on the EXACT (brute) marginals: isolates the combination step
p_anal_brute <- analytic_acat_p(matrix(p_mid, 1))

## ---- report ------------------------------------------------------------------
cat("=== Fast draws-free Cauchy vs brute-force permutation (44 blocks, N=2200) ===\n\n")
cat("[marginals] SPA two-sided p_j vs brute-force two-sided mid-p:\n")
mtab <- data.frame(rep = reps, spa = round(p_spa, 4), brute_midp = round(p_mid, 4),
                   abs_err = round(abs(p_spa - p_mid), 4))
print(mtab, row.names = FALSE)
cat(sprintf("  max |SPA - brute| over representations = %.2e\n", max(abs(p_spa - p_mid))))

cat("\n[combined p-value] three ways:\n")
cat(sprintf("  brute-force permutation Cauchy (riposte reference) = %.4f\n", cau$p.value))
cat(sprintf("  analytic ACAT on brute marginals (combine only)    = %.4f\n", p_anal_brute))
cat(sprintf("  FAST: analytic ACAT on SPA marginals (no draws)    = %.4f\n", p_fast))
cat(sprintf("  fast vs reference |diff|                           = %.4f\n",
            abs(p_fast - cau$p.value)))

cat("\n[timing]\n")
cat(sprintf("  fast path (6 SPA solves + combine) = %.4fs\n", tfast))
cat(sprintf("  brute force (B=%d draws)        = %.2fs\n", B, tbrute))
cat(sprintf("  speedup                            = %.0fx\n", tbrute / tfast))
