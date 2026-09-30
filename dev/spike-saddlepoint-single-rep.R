## spike-saddlepoint-single-rep.R --- prototype: Robinson (1982) saddlepoint tail
## for ONE riposte representation, validated against brute-force permutation.
##
## THE POINT. riposte gets its combined
## p-value by re-randomizing treatment within blocks B times and recomputing the
## statistic; the draw loop is ~72% of the cost and grows linearly in N and in B.
## For a SINGLE within-block linear statistic T = sum_i z_i s_i, the exact
## permutation tail probability is available in closed form via a saddlepoint
## approximation to the EXACT permutation cumulant generating function (CGF) --
## no draws at all, cost independent of B. This script builds that for the `raw`
## representation on the 44-block / N=2,200 example and checks it against brute
## force. The agreement check is the whole point: the SPA must track the exact
## permutation p-value to ~3 tail digits, or it is not a substitute for it.
##
## A CORRECTION TO THE BRIEF, made concrete here. The brief says riposte already
## holds the saddlepoint inputs in `riposte_sw_moments`. That returns only the
## first two cumulants (mu, Sigma). A tail built from those two alone is the
## Gaussian/MVN approximation the proposal argues against. The saddlepoint's edge
## comes from the FULL CGF -- all cumulants -- which we build below from the same
## centred scores and treated counts. Gate A checks that this CGF's first two
## derivatives at theta = 0 reproduce riposte_sw_moments exactly, so the new
## object is consistent with the old one and adds the higher cumulants on top.
##
## ON pnorm/dnorm BELOW. The spec forbids manufacturing a p-value with a probit
## (a normal approximation to the permutation distribution). The standard normal
## cdf/pdf appear here only INSIDE the Lugannani-Rice formula, which is an
## asymptotic expansion of the EXACT permutation CGF. They are special functions
## in that expansion, not a normal model for T. The distinction is the proposal's
## whole point and is preserved.

suppressMessages(devtools::load_all(".", quiet = TRUE))

## ---- 0. a seeded, balanced example: 44 blocks of 50, 25 treated each ---------
## Balanced only to keep the example clean; nothing below assumes equal blocks.
make_example <- function(n_blocks = 44L, block_size = 50L, m = 25L,
                         tau = 0.08, seed = 20260625L) {
  set.seed(seed)
  n <- n_blocks * block_size
  block <- factor(rep(seq_len(n_blocks), each = block_size))
  ## one fixed observed assignment: m treated per block, the rest control
  z <- numeric(n)
  for (ix in split(seq_len(n), block)) z[ix[seq_len(m)]] <- 1
  ## outcome: within-block noise plus a modest additive effect on the treated,
  ## so the observed statistic sits in the tail but stays brute-force estimable
  y <- rnorm(n) + tau * z
  list(y = y, z = z, block = block)
}

ex <- make_example()
## the `raw` representation, centred within block -> a single score column
S <- riposte_score_matrix(ex$y, ex$block,
                          representations = list(raw = function(y) y))
s <- S$scores[, 1]                       # length-N centred scores
blk <- ex$block

## ---- 1. the exact within-block permutation CGF for one linear statistic ------
## Per block b with centred scores s_b and treated count m_b, the statistic
## contribution sum_{i in treated} s_i is a sum over a uniform random size-m_b
## subset. Its MGF is the elementary symmetric polynomial of the tilted weights:
##   M_b(theta) = e_{m_b}(exp(theta s_b)) / C(n_b, m_b).
## Blocks are independent, so K(theta) = sum_b [ log e_{m_b}(exp(theta s_b))
##   - lchoose(n_b, m_b) ]. We compute log e_{m_b} in log space (the weights span
## many orders of magnitude at the saddlepoint) via the standard DP: process
## units one at a time, updating coefficients j = m_b .. 1 from the previous
## state. log e_k lives at index k+1.

logaddexp <- function(a, b) {                      # stable log(exp a + exp b)
  hi <- pmax(a, b); lo <- pmin(a, b)
  res <- hi + log1p(exp(lo - hi))
  inf <- is.infinite(hi)                           # both -Inf -> -Inf; +Inf -> +Inf
  res[inf] <- hi[inf]
  res
}

## log e_m of weights exp(x), x a vector of theta*s for one block
log_esym <- function(x, m) {
  logp <- c(0, rep(-Inf, m))                       # logp[k+1] = log e_k; e_0 = 1
  for (xi in x) {                                  # add one unit's weight
    ## vectorized update over k = m..1 (high-to-low so logp[k] uses old logp[k-1])
    logp[2:(m + 1)] <- logaddexp(logp[2:(m + 1)], xi + logp[1:m])
  }
  logp[m + 1]
}

## precompute per-block index sets and treated counts once (B-independent)
idx_list <- split(seq_along(s), blk)
mb_vec   <- vapply(idx_list, function(ix) sum(ex$z[ix]), numeric(1))
lognorm  <- sum(mapply(function(ix, m) lchoose(length(ix), m), idx_list, mb_vec))

cgf <- function(theta) {
  acc <- 0
  for (b in seq_along(idx_list)) {
    ix <- idx_list[[b]]; m <- mb_vec[b]
    if (m == 0L || m == length(ix)) next           # degenerate block: T_b is fixed (=0 here)
    acc <- acc + log_esym(theta * s[ix], m)
  }
  acc - lognorm
}

## derivatives by central differences; K is smooth and O(1), so a moderate h is
## accurate to far better than the 3 tail digits we need (Gate A confirms it)
.h <- 1e-3
cgf_d1 <- function(theta) (cgf(theta + .h) - cgf(theta - .h)) / (2 * .h)
cgf_d2 <- function(theta) (cgf(theta + .h) - 2 * cgf(theta) + cgf(theta - .h)) / .h^2

## ---- 2. saddlepoint solve + Lugannani-Rice upper tail ------------------------
## Solve K'(theta_hat) = t (K is convex, so K' is increasing and the root is
## unique for t strictly inside the support), then invert.
solve_saddle <- function(t) {
  mu0 <- cgf_d1(0)                                  # = permutation mean
  if (abs(t - mu0) < 1e-9) return(0)
  ## bracket: t > mu -> theta_hat > 0; expand the far end until K' overshoots t
  if (t > mu0) { lo <- 0; hi <- 1; while (cgf_d1(hi) < t) hi <- hi * 2 }
  else         { hi <- 0; lo <- -1; while (cgf_d1(lo) > t) lo <- lo * 2 }
  uniroot(function(th) cgf_d1(th) - t, c(lo, hi), tol = 1e-10)$root
}

## P(T >= t) by Lugannani-Rice; near theta_hat = 0 fall back to the (removable)
## limit, which is the normal-approximation value at the mean
spa_upper <- function(t) {
  th <- solve_saddle(t)
  if (abs(th) < 1e-6) {
    sd <- sqrt(cgf_d2(0)); return(1 - pnorm((t - cgf_d1(0)) / sd))
  }
  Kth <- cgf(th); K2 <- cgf_d2(th)
  w <- sign(th) * sqrt(2 * (th * t - Kth))
  u <- th * sqrt(K2)
  1 - pnorm(w) + dnorm(w) * (1 / u - 1 / w)
}

## ---- 3. brute-force reference (the exact target) -----------------------------
B <- 100000L
t0 <- system.time({
  draws <- riposte_block_draws(ex$z, blk, nresample = B)   # the bottleneck loop
  Tnull <- as.numeric(crossprod(draws, s))                 # T at observed (1) + draws
})[["elapsed"]]
t_obs <- Tnull[1]
emp_tail <- function(t) mean(Tnull[-1] >= t)                # population upper tail
p_perm  <- riposte_perm_pvalue(Tnull)                       # (1 + #>=)/(1 + B), observed incl.
mc_se   <- function(p) sqrt(p * (1 - p) / B)

## ---- 4. validation -----------------------------------------------------------
cat("=== Saddlepoint vs brute force: single `raw` representation ===\n")
cat(sprintf("design: %d blocks, N=%d, %d treated/block; brute B=%d (%.1fs)\n",
            nlevels(blk), length(s), unique(mb_vec)[1], B, t0))

## Gate A: the CGF's first two cumulants must match riposte_sw_moments exactly
mom <- riposte_sw_moments(matrix(s, ncol = 1), ex$z, blk)
cat("\n[Gate A] CGF cumulants vs riposte_sw_moments (closed form):\n")
cat(sprintf("  K'(0) = %.6f   mu        = %.6f   |diff| = %.2e\n",
            cgf_d1(0), mom$mu[1], abs(cgf_d1(0) - mom$mu[1])))
cat(sprintf("  K''(0)= %.6f   Sigma[1,1]= %.6f   |rel|  = %.2e\n",
            cgf_d2(0), mom$Sigma[1, 1], abs(cgf_d2(0) - mom$Sigma[1, 1]) / mom$Sigma[1, 1]))

## Gate B: the observed one-sided p-value
ts <- system.time(p_spa <- spa_upper(t_obs))[["elapsed"]]
cat("\n[Gate B] observed upper-tail p-value at t_obs =", sprintf("%.4f", t_obs), ":\n")
cat(sprintf("  brute force (pop tail)  = %.5f  (MC SE %.1e)\n", emp_tail(t_obs), mc_se(emp_tail(t_obs))))
cat(sprintf("  brute force (1+#)/(1+B) = %.5f\n", p_perm))
cat(sprintf("  saddlepoint             = %.5f  (%.4fs)\n", p_spa, ts))
cat(sprintf("  |SPA - brute(pop)|      = %.2e\n", abs(p_spa - emp_tail(t_obs))))

## Gate C: agreement across the whole distribution, including the tail. Compare
## SPA tail to empirical tail at the brute-force quantiles; report the worst
## absolute error overall and the worst RELATIVE error in the tail (p < 0.05),
## where accuracy matters and the normal approximation is known to drift.
qs <- quantile(Tnull[-1], probs = c(.5, .75, .9, .95, .975, .99, .995, .999))
tab <- data.frame(
  threshold = round(as.numeric(qs), 3),
  emp       = vapply(qs, emp_tail, numeric(1)),
  spa       = vapply(qs, spa_upper, numeric(1))
)
tab$abs_err <- abs(tab$spa - tab$emp)
tab$rel_err <- tab$abs_err / tab$emp
cat("\n[Gate C] SPA tail vs empirical tail across quantiles:\n")
print(tab, row.names = FALSE, digits = 4)
tailrows <- tab$emp < 0.05
cat(sprintf("\n  max |abs err| overall      = %.2e\n", max(tab$abs_err)))
cat(sprintf("  max |rel err| in tail(<.05) = %.2e\n", max(tab$rel_err[tailrows])))

## ---- 5. how the cost scales: SPA is free of B and cheap in N ------------------
cat("\n[scaling] saddlepoint cost is independent of B and grows only with N:\n")
for (nb in c(44L, 440L, 4400L)) {
  ex2 <- make_example(n_blocks = nb)
  S2  <- riposte_score_matrix(ex2$y, ex2$block, representations = list(raw = function(y) y))
  s2  <- S2$scores[, 1]; blk2 <- ex2$block
  idx2 <- split(seq_along(s2), blk2); mb2 <- vapply(idx2, function(ix) sum(ex2$z[ix]), numeric(1))
  ln2  <- sum(mapply(function(ix, m) lchoose(length(ix), m), idx2, mb2))
  ## rebind the closures' data via a local CGF (kept inline to avoid global state)
  cgf2 <- function(theta) { acc <- 0
    for (b in seq_along(idx2)) { ix <- idx2[[b]]; m <- mb2[b]
      if (m == 0L || m == length(ix)) next; acc <- acc + log_esym(theta * s2[ix], m) }
    acc - ln2 }
  d1 <- function(th) (cgf2(th + .h) - cgf2(th - .h)) / (2 * .h)
  d2 <- function(th) (cgf2(th + .h) - 2 * cgf2(th) + cgf2(th - .h)) / .h^2
  ss <- function(t) { if (t > 0) { hi <- 1; while (d1(hi) < t) hi <- hi * 2; lo <- 0 }
    else { lo <- -1; while (d1(lo) > t) lo <- lo * 2; hi <- 0 }
    uniroot(function(th) d1(th) - t, c(lo, hi), tol = 1e-10)$root }
  tt <- as.numeric(crossprod(ex2$z, s2))            # observed statistic
  el <- system.time({ th <- ss(tt); w <- sign(th) * sqrt(2 * (th * tt - cgf2(th)))
    u <- th * sqrt(d2(th)); p <- 1 - pnorm(w) + dnorm(w) * (1 / u - 1 / w) })[["elapsed"]]
  cat(sprintf("  N=%6d (%4d blocks): SPA p=%.4f in %.3fs\n", length(s2), nb, p, el))
}
cat("\n(Pure-R DP; a compiled inner loop makes each SPA call sub-millisecond.\n",
    " The algorithmic point is the B-free O(sum_b n_b m_b) cost.)\n", sep = "")
