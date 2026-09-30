## stephenson.R --- Stephenson rank scores as outcome representations.
##
## THE POINT. The Stephenson rank test targets the UPPER TAIL: with score
## choose(rank - 1, r - 1) for a unit at ascending within-block rank, a larger
## tuning value r upweights the top ranks, so the statistic responds to effects
## concentrated in a few units helped a lot (the "fraction helped by at least c"
## alternative) rather than a uniform shift. r = 2 is the Wilcoxon score.
##
## A single r commits to one tail weighting. The COMBINED Stephenson test hedges
## across several r. In riposte that is automatic: each r is just another
## within-block representation, so riposte_stephenson_reps() returns a set of
## representations that the existing quadratic/Cauchy/screen combine, referred to
## the randomization distribution with the exact closed-form moments. This is the
## combined Stephenson test of Kim, Su, Bowers, and Li, expressed in riposte's
## representation framework rather than duplicated from CMRSS.

#' Stephenson rank scores
#'
#' For an outcome vector, the Stephenson score of a unit at ascending rank `j` is
#' `choose(j - 1, r - 1)`. Ties take average ranks. Larger `r` concentrates weight
#' on the top ranks; `r=2` reduces to the (rank-1) Wilcoxon score.
#'
#' @param y numeric vector (a single block's outcomes).
#' @param r the Stephenson tuning value, a whole number; `r=2` is the Wilcoxon
#'   (rank - 1) score and `r=1` gives a constant (uninformative) score. Must be
#'   a whole number: a fractional `r` would silently round inside `choose()`.
#' @return numeric vector of scores the same length as `y`.
#' @export
riposte_stephenson_scores <- function(y, r) {
  if (length(r) != 1L || !is.finite(r) || r < 1 || r != round(r))
    stop("`r` must be a single whole number >= 1.", call. = FALSE)
  ranks <- rank(y, ties.method = "average")
  choose(ranks - 1, r - 1)
}

#' Stephenson rank-score representations at several tuning values
#'
#' Returns a named list of representation functions, one per tuning value `r`,
#' suitable for [riposte_test()] or [riposte_components()]. Combining them (the
#' default quadratic or the screen) gives the combined Stephenson rank test, which
#' targets the upper tail across a range of weightings instead of committing to
#' one. Add the default representations with `c(riposte_reps_default(),
#' riposte_stephenson_reps())` to combine distribution-shape and upper-tail views
#' together.
#'
#' @param r vector of whole-number tuning values; `c(2, 6, 10)` by default
#'   (Wilcoxon plus two heavier upper-tail weightings), following the proposal.
#' @return a named list of representation functions, named `S2`, `S6`, ...
#' @seealso [riposte_stephenson_scores()], [riposte_reps_default()].
#' @export
riposte_stephenson_reps <- function(r = c(2, 6, 10)) {
  stats::setNames(
    lapply(r, function(rr) {
      force(rr)
      function(y) riposte_stephenson_scores(y, rr)
    }),
    paste0("S", r))
}
