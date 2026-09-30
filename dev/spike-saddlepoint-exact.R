## spike-saddlepoint-exact.R --- validate the single-representation saddlepoint
## against an EXACT permutation tail (full within-block enumeration, zero Monte
## Carlo noise) and exercise the log-space tail in the deep tail.
##
## THE POINT. The brute-force spike could only confirm agreement to within Monte
## Carlo error (the reference's own SE swamps the comparison past p ~ 0.01). Here
## the blocks are small enough to enumerate the entire permutation orbit, so the
## permutation tail is EXACT. Any gap is then pure saddlepoint error, and we can
## probe it down to ~1/orbit ~ 1e-5 -- where it actually matters and where the
## normal approximation is known to drift.

source("dev/saddlepoint-core.R")
suppressMessages(devtools::load_all(".", quiet = TRUE))

## ---- small enumerable design: 4 blocks of 6, 3 treated each ------------------
## Orbit size = C(6,3)^4 = 20^4 = 160,000, so the exact tail resolves to 6.25e-6.
set.seed(20260625L)
n_blocks <- 4L; block_size <- 6L; m <- 3L
N <- n_blocks * block_size
block <- factor(rep(seq_len(n_blocks), each = block_size))
z <- numeric(N); for (ix in split(seq_len(N), block)) z[ix[seq_len(m)]] <- 1
y <- rnorm(N) + 0.6 * z

S <- riposte_score_matrix(y, block, representations = list(raw = function(yy) yy))
s <- S$scores[, 1]
spa <- make_spa(s, z, block)

## ---- exact permutation distribution of T by convolution ----------------------
## T = sum_b T_b; within block b, T_b ranges over all C(n_b, m_b) subset sums,
## each equally likely. The orbit is the product across blocks, so the exact
## distribution of T is the convolution of the per-block subset-sum vectors.
idx_list <- split(seq_along(s), block)
per_block <- lapply(idx_list, function(ix) combn(s[ix], m, FUN = sum))   # all subset sums
T_all <- Reduce(function(acc, v) as.vector(outer(acc, v, "+")), per_block)
exact_tail <- function(t) mean(T_all >= t)                              # EXACT, no MC
## mid-tail: the discrete tail is P(T>=t); a continuous approximation targets the
## SMOOTHED tail, ~ halfway between P(T>=t) and P(T>t). Comparing to the mid value
## removes the half-atom discreteness bias that otherwise inflates the deep-tail
## relative error against a continuous formula.
exact_mid  <- function(t) (mean(T_all >= t) + mean(T_all > t)) / 2
cat(sprintf("orbit size = %d (exact tail resolves to %.2e); support [%.2f, %.2f]\n",
            length(T_all), 1 / length(T_all), spa$supmin, spa$supmax))

## ---- Gate A: CGF cumulants vs closed-form moments ----------------------------
mom <- riposte_sw_moments(matrix(s, ncol = 1), z, block)
cat(sprintf("\n[Gate A] K'(0)=%.6f vs mu=%.6f (|d| %.1e);  K''(0)=%.6f vs Sigma=%.6f (|rel| %.1e)\n",
            spa$mu, mom$mu[1], abs(spa$mu - mom$mu[1]),
            spa$sigma2, mom$Sigma[1, 1], abs(spa$sigma2 - mom$Sigma[1, 1]) / mom$Sigma[1, 1]))

## ---- Gate C': SPA vs EXACT tail across the distribution, into the deep tail ---
## thresholds chosen so the exact tail hits target levels (type=1: land on actual
## support points, so exact_tail returns those levels cleanly)
targets <- c(.1, .05, .01, .005, .001, .0005, .0001, .00005)
thr <- as.numeric(quantile(T_all, probs = 1 - targets, type = 1))
tab <- data.frame(
  target_p  = targets,
  threshold = round(thr, 3),
  exact     = vapply(thr, exact_tail, numeric(1)),
  mid       = vapply(thr, exact_mid,  numeric(1)),
  spa       = vapply(thr, spa$upper, numeric(1))
)
tab$abs_err     <- abs(tab$spa - tab$exact)
tab$rel_err     <- tab$abs_err / tab$exact
tab$rel_err_mid <- abs(tab$spa - tab$mid) / tab$mid    # vs the de-discretised tail
cat("\n[Gate C'] SPA tail vs EXACT enumerated tail (no Monte Carlo noise):\n")
print(tab, row.names = FALSE, digits = 4)
cat(sprintf("\n  max |abs err| overall                         = %.2e\n", max(tab$abs_err)))
cat(sprintf("  max |rel err| vs discrete tail (deepest p~5e-5)= %.2e\n", max(tab$rel_err)))
cat(sprintf("  max |rel err| vs mid (half-atom bias removed)  = %.2e\n", max(tab$rel_err_mid)))
ld <- vapply(thr, spa$upper_log, numeric(1))
cat(sprintf("  direct vs log-space tail agree to             = %.1e\n",
            max(abs(tab$spa - exp(ld)))))

## ---- interior deep-tail: log p where simulation cannot reach ------------------
## A threshold just below the support maximum has a tiny but NONZERO exact tail
## (a real, possible event), far smaller than 1/B for any feasible brute-force B.
## The SPA returns it directly in log space.
thr_deep <- sort(unique(T_all), decreasing = TRUE)[3]   # 3rd-largest attainable value
cat(sprintf("\n[interior deep tail] at t = %.3f (just below supmax = %.3f):\n",
            thr_deep, spa$supmax))
cat(sprintf("  exact enumerated tail     = %.3e  (only %d of %d orbit elements)\n",
            exact_tail(thr_deep), sum(T_all >= thr_deep), length(T_all)))
cat(sprintf("  saddlepoint  P(T>=t)      = %.3e  (log10 = %.2f)\n",
            exp(spa$upper_log(thr_deep)), spa$upper_log(thr_deep) / log(10)))
