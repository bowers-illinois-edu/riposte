## saddlepoint-core.R --- reusable Robinson (1982) saddlepoint tail for ONE
## within-block linear statistic. Sourced by the dev spike scripts; not package
## code yet. `make_spa(s, z, block)` returns the CGF and tail functions for the
## linear statistic T = sum_i z_i s_i referred to the within-block permutation
## null (each block's treated count fixed). The construction is exact in the
## sense that the CGF is the EXACT permutation CGF (all cumulants); only the tail
## inversion (Lugannani-Rice) is an asymptotic approximation to that exact CGF.

## stable log(exp a + exp b), vectorized
logaddexp <- function(a, b) {
  hi <- pmax(a, b); lo <- pmin(a, b)
  res <- hi + log1p(exp(lo - hi))
  inf <- is.infinite(hi)                       # both -Inf -> -Inf; +Inf -> +Inf
  res[inf] <- hi[inf]
  res
}

## log of the elementary symmetric polynomial e_m of weights exp(x), computed in
## log space because the weights span many orders of magnitude at the saddlepoint.
## DP: process units one at a time; log e_k lives at index k+1; update high-to-low.
log_esym <- function(x, m) {
  logp <- c(0, rep(-Inf, m))                   # log e_0 = 0
  for (xi in x) logp[2:(m + 1)] <- logaddexp(logp[2:(m + 1)], xi + logp[1:m])
  logp[m + 1]
}

## Build the saddlepoint object for the linear statistic with centred scores `s`,
## observed assignment `z`, and `block`. Returns a list of closures.
make_spa <- function(s, z, block) {
  s <- as.numeric(s); z <- as.numeric(z); block <- as.factor(block)
  idx_list <- split(seq_along(s), block)
  mb  <- vapply(idx_list, function(ix) sum(z[ix]), numeric(1))
  nb  <- lengths(idx_list)

  ## exact support of T: per block the largest/smallest the contribution can be
  ## is the sum of the top/bottom m_b scores. T can never exceed supmax or fall
  ## below supmin, so a tail asked outside the support is 0 or 1 with no
  ## saddlepoint (the equation K'(theta)=t has no finite root there).
  supmax <- sum(mapply(function(ix, m) sum(sort(s[ix], decreasing = TRUE)[seq_len(m)]),
                       idx_list, mb))
  supmin <- sum(mapply(function(ix, m) sum(sort(s[ix])[seq_len(m)]), idx_list, mb))

  ## K(theta) = sum_b log M_b(theta). A block with m_b in {0, n_b} contributes a
  ## deterministic amount (0 if scores are within-block centred); the rest use the
  ## elementary-symmetric MGF normalised by C(n_b, m_b).
  cgf <- function(theta) {
    acc <- 0
    for (b in seq_along(idx_list)) {
      ix <- idx_list[[b]]; m <- mb[b]; n <- nb[b]
      if (m == 0L) next
      if (m == n)  { acc <- acc + theta * sum(s[ix]); next }   # deterministic block
      acc <- acc + log_esym(theta * s[ix], m) - lchoose(n, m)
    }
    acc
  }

  ## central differences: K is smooth and O(1), so a moderate step is accurate to
  ## far better than the tail digits we need (checked against riposte_sw_moments)
  hstep <- 1e-3
  d1 <- function(th) (cgf(th + hstep) - cgf(th - hstep)) / (2 * hstep)
  d2 <- function(th) (cgf(th + hstep) - 2 * cgf(th) + cgf(th - hstep)) / hstep^2

  ## solve K'(theta_hat) = t (K convex => K' increasing => unique interior root)
  mu0 <- d1(0)
  ## solve K'(theta_hat) = t for t strictly inside the support; cap the bracket so
  ## a t at the support edge cannot send hi to infinity (NaN). Callers short-
  ## circuit t outside the support before reaching here.
  solve_saddle <- function(t) {
    if (abs(t - mu0) < 1e-9) return(0)
    if (t > mu0) { lo <- 0; hi <- 1; while (d1(hi) < t && hi < 1e3) hi <- hi * 2 }
    else         { hi <- 0; lo <- -1; while (d1(lo) > t && lo > -1e3) lo <- lo * 2 }
    uniroot(function(th) d1(th) - t, c(lo, hi), tol = 1e-10)$root
  }

  ## Lugannani-Rice upper tail P(T >= t), direct form (underflows in the deep tail)
  upper <- function(t) {
    if (t > supmax) return(0); if (t <= supmin) return(1)
    th <- solve_saddle(t)
    if (abs(th) < 1e-6) return(stats::pnorm((t - mu0) / sqrt(d2(0)), lower.tail = FALSE))
    w <- sign(th) * sqrt(2 * (th * t - cgf(th))); u <- th * sqrt(d2(th))
    1 - stats::pnorm(w) + stats::dnorm(w) * (1 / u - 1 / w)
  }

  ## log P(T >= t), stable in the deep tail. Factor out dnorm(w):
  ##   tail = phi(w) * [ R(w) + 1/u - 1/w ],  R(w) = Phibar(w)/phi(w) the Mills ratio.
  ## R(w) is computed as exp(log Phibar - log phi), so it never underflows. Intended
  ## for the upper tail (t >= mu0, th >= 0).
  upper_log <- function(t) {
    if (t > supmax) return(-Inf); if (t <= supmin) return(0)
    th <- solve_saddle(t)
    if (abs(th) < 1e-6)
      return(stats::pnorm((t - mu0) / sqrt(d2(0)), lower.tail = FALSE, log.p = TRUE))
    w <- sign(th) * sqrt(2 * (th * t - cgf(th))); u <- th * sqrt(d2(th))
    R <- exp(stats::pnorm(w, lower.tail = FALSE, log.p = TRUE) - stats::dnorm(w, log = TRUE))
    stats::dnorm(w, log = TRUE) + log(R + 1 / u - 1 / w)
  }

  ## two-sided tail about the permutation mean, P(|T - mu0| >= |t - mu0|), the
  ## continuous-SPA analogue of riposte's two-sided mid-p marginal. The upper()
  ## short-circuits handle thresholds outside the support (one side contributes 0).
  two_sided <- function(t) {
    dev <- abs(t - mu0)
    min(1, upper(mu0 + dev) + (1 - upper(mu0 - dev)))
  }

  list(cgf = cgf, d1 = d1, d2 = d2, solve_saddle = solve_saddle,
       upper = upper, upper_log = upper_log, two_sided = two_sided,
       mu = mu0, sigma2 = d2(0), supmin = supmin, supmax = supmax)
}
