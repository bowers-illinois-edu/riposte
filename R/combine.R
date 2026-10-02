## combine.R --- the three combinations.
##
## THE POINT. From the representation linear statistics, riposte forms three
## combinations, all referred to the randomization distribution:
##   quadratic  (T - mu)' Sigma^+ (T - mu) using the FULL permutation covariance
##              (closed-form mu, Sigma; Moore-Penrose inverse via MASS::ginv).
##              The STATISTIC uses the exact closed-form moments (not estimated
##              from the draws, which made the proposal's version anti-conservative
##              by O(1/nresample)); its P-VALUE is the share of the randomization
##              distribution at least as large as observed. Pools overlapping
##              evidence.
##   cauchy     mean of the pole-aware Cauchy transform of the representations'
##              two-sided MID-P permutation p-values; uses only marginal
##              calibration, ignores the covariance. Mid-p keeps the transform off
##              its poles --- no clamp, no probit.
##   max        the largest standardized representation; a third comparison.
##
## The Cauchy and the max also take a one-sided alternative. "greater" counts
## only treated-score sums above their permutation mean, "less" only sums below
## it: the Cauchy averages the transforms of one-sided mid-p values, and the max
## takes the largest standardized sum in the named direction. The quadratic
## squares every departure, so it has no one-sided form. With scores that rise
## with the outcome, "less" asks whether treated outcomes are lower (harm).
##
## All three read the SAME shared draws (riposte_perm_stats), so they are compared
## on one randomization distribution. Only the references are Monte-Carlo; the
## quadratic's moments and the max's scales are closed-form.

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- combinations on a statistics matrix ------------------------------------
## Tmat is (1 + nresample) x representations, row 1 the observed assignment. The
## unadjusted path builds it as crossprod(draws, fixed scores) and supplies the
## closed-form moments; the adjusted path builds it from per-permutation residuals
## and supplies Monte-Carlo moments. These three functions are the shared core.

riposte_quadratic_from_T <- function(Tmat, mu, Sigma) {
  Sigma_inv <- MASS::ginv(Sigma)            # handles the rank deficiency the screen guards against
  cen <- sweep(Tmat, 2, mu)
  q <- rowSums((cen %*% Sigma_inv) * cen)
  list(statistic = unname(q[1]), p.value = riposte_perm_pvalue(q), df = qr(Sigma)$rank)
}

riposte_cauchy_from_T <- function(Tmat, alternative = "two.sided") {
  ## mid-p of each rep at each assignment, in the requested direction
  P <- apply(Tmat, 2, riposte_midp, alternative = alternative)
  Tc <- rowMeans(riposte_acat_term(P))
  list(statistic = unname(Tc[1]), p.value = riposte_perm_pvalue(Tc),
       component_midp = stats::setNames(P[1, ], colnames(Tmat)))
}

riposte_max_from_T <- function(Tmat, mu, sdv, alternative = "two.sided") {
  sdv[sdv == 0] <- 1
  Z <- sweep(sweep(Tmat, 2, mu), 2, sdv, "/")
  ## orient so that a large value is evidence in the requested direction; for
  ## "less" the statistic is the largest standardized shortfall
  Zdir <- switch(alternative, two.sided = abs(Z), greater = Z, less = -Z)
  mx <- apply(Zdir, 1, max)
  list(statistic = unname(mx[1]), p.value = riposte_perm_pvalue(mx),
       which = colnames(Tmat)[which.max(Zdir[1, ])])
}

#' Quadratic (energy/distance) omnibus
#'
#' @param scores within-block-centred score matrix.
#' @param z 0/1 treatment vector.
#' @param block block factor.
#' @param moments optional precomputed [riposte_sw_moments()] result (shared).
#' @param draws optional precomputed [riposte_block_draws()] matrix (shared).
#' @param nresample re-randomizations, used only when `draws` is `NULL`.
#' @return a list with the observed `statistic`, its permutation `p.value`, and
#'   `df` (the rank of Sigma, the degrees of freedom of the quadratic form).
#' @keywords internal
#' @noRd
riposte_quadratic <- function(scores, z, block, moments = NULL, draws = NULL,
                              nresample = 1999L) {
  scores <- as.matrix(scores)
  if (ncol(scores) == 0L)
    return(list(statistic = NA_real_, p.value = NA_real_, df = 0L))
  moments <- moments %||% riposte_sw_moments(scores, z, block)
  Tmat <- riposte_perm_stats(scores, z, block, nresample, draws)$stats
  riposte_quadratic_from_T(Tmat, moments$mu, moments$Sigma)
}

#' Cauchy (ACAT) combination of the representations' mid-p permutation p-values
#'
#' @inheritParams riposte_quadratic
#' @param alternative `"two.sided"`, `"greater"`, or `"less"`: the direction of
#'   the mid-p values that are combined.
#' @return a list with the observed `statistic`, its permutation `p.value`, and
#'   `component_midp` (each representation's observed mid-p value).
#' @keywords internal
#' @noRd
riposte_cauchy <- function(scores, z, block, draws = NULL, nresample = 1999L,
                           alternative = "two.sided") {
  scores <- as.matrix(scores)
  if (ncol(scores) == 0L)
    return(list(statistic = NA_real_, p.value = NA_real_,
                component_midp = numeric(0)))
  Tmat <- riposte_perm_stats(scores, z, block, nresample, draws)$stats
  riposte_cauchy_from_T(Tmat, alternative)
}

#' Max combination (largest standardized representation)
#'
#' @inheritParams riposte_quadratic
#' @param alternative `"two.sided"` (largest absolute standardized sum),
#'   `"greater"` (largest), or `"less"` (largest in the negative direction).
#' @return a list with the observed `statistic`, its permutation `p.value`, and
#'   `which` (the representation attaining the observed maximum).
#' @keywords internal
#' @noRd
riposte_max <- function(scores, z, block, moments = NULL, draws = NULL,
                        nresample = 1999L, alternative = "two.sided") {
  scores <- as.matrix(scores)
  if (ncol(scores) == 0L)
    return(list(statistic = NA_real_, p.value = NA_real_, which = NA_character_))
  moments <- moments %||% riposte_sw_moments(scores, z, block)
  Tmat <- riposte_perm_stats(scores, z, block, nresample, draws)$stats
  riposte_max_from_T(Tmat, moments$mu, sqrt(diag(moments$Sigma)), alternative)
}
