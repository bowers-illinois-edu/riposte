## poly-ranks.R --- polynomial rank scores as outcome representations.
##
## THE POINT. Like the Stephenson score (stephenson.R), the polynomial rank score
## puts its weight on the top of each block's ranking, so its statistic responds
## to a treatment that raises the outcomes of a few units a great deal. The unit
## at within-block rank k of n_b gets (k / (n_b + 1))^(zeta - 1). The paper (Bowers
## and Burton) prefers this form to Stephenson's choose(k - 1, zeta - 1) because
## it stays in [0, 1) in every block: the binomial coefficient reaches about 3.9e13
## at rank 50 when zeta = 22, so blocks of different sizes would enter on scales
## that differ by orders of magnitude, and in the paper's code that spread broke
## the quadratic combination's covariance inverse.
##
## Each zeta is one more representation, so the existing combinations (quadratic,
## Cauchy, screen, max) combine several of them with no new machinery. The
## default zetas reproduce the paper's rank-score tables; the vignette shows how
## to choose others for a given block size.
##
## tail = "lower" weights the BOTTOM of each block instead, for a treatment that
## lowers the outcomes of a few units a great deal (harm). It ranks from the top,
## (n_b + 1 - k) / (n_b + 1), and takes the negative, so that the score still RISES
## with the outcome. Every built-in monotone score (raw, rank, both tails) then
## rises with the outcome, and alternative = "less" in riposte_test() means
## "treated outcomes lower" for all of them at once.

#' Polynomial rank scores
#'
#' For an outcome vector of length `n`, the polynomial rank score of a unit at
#' ascending rank `k` is `(k / (n + 1))^(zeta - 1)`, with ties at their average
#' rank. This is the score of Kim, Su, Bowers, and Li (arXiv 2605.08027). At
#' `zeta = 2` it is the rank divided by `n + 1`, which gives the Wilcoxon rank-sum
#' test. Larger `zeta` puts more of the block's total score on its top ranks:
#' treating `k / (n + 1)` as uniform on (0, 1), the score-weighted average
#' distance from the top of the block is about `n / (zeta + 1)` units.
#'
#' Unlike [riposte_stephenson_scores()], the score lies in `[0, 1)` whatever the
#' block size, so blocks of different sizes enter on the same scale.
#'
#' With `tail = "lower"` the score weights the bottom of the block instead: the
#' unit at ascending rank `k` gets `-((n + 1 - k) / (n + 1))^(zeta - 1)`, which
#' lies in `(-1, 0]`. The minus sign makes the score rise with the outcome, as
#' the upper-tail score does, so that with `alternative = "less"` in
#' [riposte_test()] both tails ask whether treated outcomes are lower. Without
#' ties, the lower-tail score is `-riposte_poly_scores(-y, zeta)`.
#'
#' @param y numeric vector (a single block's outcomes).
#' @param zeta the tuning value, a single number `>= 1`. It need not be a whole
#'   number. `zeta = 1` gives a constant, which carries no information.
#' @param tail `"upper"` (default) weights the top of the block; `"lower"`
#'   weights the bottom.
#' @return numeric vector of scores the same length as `y`.
#' @seealso [riposte_poly_reps()], [riposte_stephenson_scores()].
#' @export
riposte_poly_scores <- function(y, zeta, tail = c("upper", "lower")) {
  tail <- match.arg(tail)
  if (length(zeta) != 1L || !is.finite(zeta) || zeta < 1)
    stop("`zeta` must be a single finite number >= 1.", call. = FALSE)
  ranks <- rank(y, ties.method = "average")
  n1 <- length(y) + 1
  if (tail == "upper") (ranks / n1)^(zeta - 1)
  else -((n1 - ranks) / n1)^(zeta - 1)
}

#' Polynomial rank-score representations at several tuning values
#'
#' Returns a named list of representation functions, one per value of `zeta`,
#' for [riposte_test()] or [riposte_components()]. Combining them gives a test
#' that responds to gains concentrated at the top of each block across several
#' weightings instead of one.
#'
#' The default, `c(2, 7, 12, 17, 22)`, is the set in Bowers and Burton's
#' rank-score tables, chosen for blocks of 50 units with 25 treated. For other
#' designs, supply your own values; the vignette shows how to choose them. In
#' brief: the largest `zeta` near `n_b / m - 1` targets gains confined to about
#' `m` units per block; and two scores at `zeta` and `zeta'` correlate at about
#' `2 sqrt(u v) / (u + v)` with `u = 2 zeta - 1` and `v = 2 zeta' - 1`, so
#' values whose `2 zeta - 1` grow by a constant factor (for example `zeta = 2, 5,
#' 14, 41`, where `2 zeta - 1` triples) are evenly spaced, while equal steps in
#' `zeta` give neighbours whose correlation approaches 1.
#'
#' `zeta = 2` is the rank rescaled. When it is added to
#' [riposte_reps_default()], which already holds the rank, leave it out:
#' `c(riposte_reps_default(), riposte_poly_reps(c(7, 12, 17, 22)))`.
#'
#' With `tail = "lower"` the scores weight the bottom of each block, for a
#' treatment that lowers the outcomes of a few units; see
#' [riposte_poly_scores()] and `vignette("harm")`.
#'
#' @param zeta vector of tuning values, each `>= 1`.
#' @param tail `"upper"` (default) or `"lower"`, passed to
#'   [riposte_poly_scores()].
#' @return a named list of representation functions, named `poly2`, `poly7`, ...
#'   for the upper tail and `polylow2`, `polylow7`, ... for the lower tail.
#' @seealso [riposte_poly_scores()], [riposte_stephenson_reps()].
#' @export
riposte_poly_reps <- function(zeta = c(2, 7, 12, 17, 22), tail = c("upper", "lower")) {
  tail <- match.arg(tail)
  ## check every value now, not when the first block is scored
  for (zz in zeta) riposte_poly_scores(1, zz)
  stats::setNames(
    lapply(zeta, function(zz) {
      force(zz)
      function(y) riposte_poly_scores(y, zz, tail = tail)
    }),
    paste0(if (tail == "upper") "poly" else "polylow", zeta))
}
