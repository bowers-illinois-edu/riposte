## perm.R --- the permutation engine (the hot loop).
##
## THE POINT. Every Monte-Carlo combination (Cauchy, max, and the quadratic's
## p-value reference) needs the same object: the linear statistics of all
## representations evaluated at the observed assignment and at many within-block
## re-randomizations. We compute it once and share it, so the combinations are
## compared on one randomization distribution.
##
## THE SEAM FOR A FUTURE C++/Rust BACKEND. This is the only part of riposte that
## loops over permutations; the quadratic and the screen are closed-form and do
## not. So this function --- draw assignments, form T = G' S --- is where a
## compiled backend would pay off, and it is deliberately the whole permutation
## surface. Any such backend must reproduce the contract tested here: column/row
## 1 is the observed assignment, and the empirical moments of the draws must
## match the closed-form Strasser-Weber moments (a distributional check, since
## RNG streams will not match across languages). Keep that test when the backend
## lands; it, not bit-for-bit equality, is what guarantees correctness.

#' Linear statistics across the observed assignment and re-randomizations
#'
#' Builds the matrix of representation linear statistics `T_{j,c} = sum_i G_{i,j}
#' s_{i,c}` for the observed assignment (row 1) and `nresample` within-block
#' re-randomizations (rows 2..). This is the shared randomization distribution the
#' Monte-Carlo combinations refer to.
#'
#' @param scores numeric matrix of within-block-centred scores (units x
#'   representations), e.g. `riposte_score_matrix()$scores`.
#' @param z 0/1 treatment vector.
#' @param block factor of block labels.
#' @param nresample number of re-randomizations.
#' @param draws optional precomputed assignment matrix from
#'   [riposte_block_draws()] (units x (1 + nresample), column 1 = observed); if
#'   supplied, `nresample` is ignored and the draws are reused (this is how the
#'   combinations share one set).
#' @return a list with `stats` (a (1 + nresample) x representations matrix, row 1
#'   observed) and `draws` (the assignment matrix used).
#' @keywords internal
#' @noRd
riposte_perm_stats <- function(scores, z, block, nresample = 1999L, draws = NULL) {
  scores <- as.matrix(scores)
  if (is.null(draws)) draws <- riposte_block_draws(z, block, nresample)
  ## T = G' S: each column of G is an assignment, each column of S a
  ## representation; crossprod gives assignments x representations
  stats <- crossprod(draws, scores)
  colnames(stats) <- colnames(scores)
  list(stats = stats, draws = draws)
}

#' Upper-tail permutation p-value with the observed value included
#'
#' For a statistic evaluated across assignments (`stat[1]` observed, the rest
#' draws), the share at least as large as the observed value. Because the
#' observed value is included in both numerator and denominator, this is the
#' exactly valid `(1 + #draws >= observed) / (1 + nresample)`.
#'
#' @param stat numeric vector of a statistic across assignments, observed first.
#' @return a single p-value in (0, 1].
#' @keywords internal
#' @noRd
riposte_perm_pvalue <- function(stat) {
  mean(stat >= stat[1])
}
