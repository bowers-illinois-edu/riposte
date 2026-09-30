## representations.R --- the outcome representations (transformations).
##
## THE POINT. riposte builds several within-block representations of an outcome
## Y and combines the evidence each carries. A representation is just a function
## that maps a block's outcome vector to a per-unit score vector of the same
## length: raw values, ranks, a unit's mean distance to the others, its max
## distance, a bounded transform, or anything the user supplies. The combined
## tests then centre each score within block and sum over the treated units (the
## Strasser-Weber linear-statistic form).
##
## Two design commitments live here:
##   1. The set of representations is OPEN. riposte_reps_default() returns the
##      six of Bowers and Burton's working paper (Huber's psi sixth; the NSF
##      proposal's raw-tanh set is riposte_reps_proposal()), but
##      riposte_test()/riposte_components() accept any
##      named list of functions, so a user can drop in log distances, a winsorised
##      mean, a kernel score, or replace the set entirely.
##   2. The default distance representations are computed in a numerically stable
##      way. The mean L1 distance is translation invariant, so we subtract the
##      median before summing; that keeps the partial sums O(spread) instead of
##      O(|location|). The payoff: the fast O(n log n) closed form then matches
##      the definitional O(n^2) sum-of-distances on the SAME input to machine
##      precision at any location, whereas the plain cumsum-then-subtract form
##      (the manytestsr C++) adds summation-cancellation error --- about 3e-4 on an
##      outcome offset by 1e12 (dollars, populations, raw test points), some 5000x
##      the unavoidable input-rounding floor. Centring cannot recover precision
##      already lost when the outcome is STORED as offset + spread in a double
##      (no algorithm can); it only avoids ADDING error of its own.

#' Mean pairwise L1 distance of each unit to the others (within a block)
#'
#' For an outcome vector `y` of length n, returns, for each unit i, the mean
#' absolute distance `mean_{j != i} |y_i - y_j|` with denominator `n - 1`. This
#' is the energy/distance representation: units far from the rest of their block
#' (in either tail) score high.
#'
#' The computation is the sorted closed form (no n-by-n distance matrix), made
#' numerically stable by first subtracting the median. L1 distances are
#' translation invariant, so the subtraction changes nothing about the answer; it
#' keeps the prefix sums small so the fast form matches the definitional
#' sum-of-distances on the same input to machine precision even when `y` carries a
#' large location (dollars, populations). It does not recover precision already
#' lost when `y` is stored as a large offset plus a small spread --- no algorithm
#' can --- but unlike a plain cumsum-then-subtract it adds no error of its own.
#'
#' @param y numeric vector (a single block's outcomes).
#' @return numeric vector the same length as `y`. A block of length < 2 has no
#'   distances and returns zeros.
#' @export
riposte_mean_dist <- function(y) {
  n <- length(y)
  if (n < 2L) return(numeric(n))
  ## translation invariant: centre to keep summation well-conditioned
  yc  <- y - stats::median(y)
  ord <- order(yc)
  ys  <- yc[ord]
  csum <- cumsum(ys)
  k         <- seq_len(n)
  below_sum <- csum - ys          # sum of strictly-smaller elements
  above_sum <- csum[n] - csum     # sum of strictly-larger elements (csum[n] ~ 0)
  below_cnt <- k - 1L
  above_cnt <- n - k
  d_sorted  <- (ys * below_cnt - below_sum) + (above_sum - ys * above_cnt)
  d <- numeric(n)
  d[ord] <- d_sorted / (n - 1)
  d
}

#' Mean pairwise L1 distance computed on the within-block ranks
#'
#' Like [riposte_mean_dist()] but applied to the mid-ranks of `y`, so it measures
#' how far each unit sits from the others in rank space. Robust to the outcome's
#' scale and to outliers.
#'
#' @inheritParams riposte_mean_dist
#' @return numeric vector the same length as `y`.
#' @export
riposte_mean_rank_dist <- function(y) {
  riposte_mean_dist(riposte_rank(y))
}

#' Maximum distance of each unit to the block extremes
#'
#' For each unit, `max(y_i - min(y), max(y) - y_i)`: how far the unit is from the
#' nearer end of its block's range. Translation invariant.
#'
#' @inheritParams riposte_mean_dist
#' @return numeric vector the same length as `y`.
#' @export
riposte_max_dist <- function(y) {
  if (length(y) < 1L) return(numeric(0))
  pmax(y - min(y), max(y) - y)
}

#' Mid-ranks of a block's outcomes
#'
#' Average ranks (ties resolved by averaging), matching the rank representation
#' used by the energy combination.
#'
#' @inheritParams riposte_mean_dist
#' @return numeric vector the same length as `y`.
#' @export
riposte_rank <- function(y) {
  rank(y, ties.method = "average")
}

#' Bounded (tanh) transform of a block's outcomes, robustly standardised
#'
#' A bounded representation: `tanh((y - median(y)) / s)`, where `s` is a robust
#' scale (the MAD, falling back to the SD and then to 1 when those are zero).
#'
#' The standardisation matters. Raw `tanh(y)` saturates to +/-1 once `|y|`
#' exceeds about 3, so on an outcome measured in dollars or raw test points every
#' value maps to +/-1 and the representation carries no information. Centring and
#' scaling first keeps the transform sensitive in the body of each block's
#' distribution while still bounding the influence of extreme values. To
#' reproduce the proposal's raw `tanh(y)` exactly, use [riposte_reps_proposal()].
#'
#' This transform is not in either built-in set: [riposte_reps_default()] uses
#' [riposte_huber()] in the sixth place. Add it to a list of representations to
#' use it.
#'
#' @inheritParams riposte_mean_dist
#' @return numeric vector the same length as `y`.
#' @export
riposte_tanh <- function(y) {
  s <- stats::mad(y)
  if (!is.finite(s) || s <= 0) s <- stats::sd(y)
  if (!is.finite(s) || s <= 0) s <- 1
  tanh((y - stats::median(y)) / s)
}

#' Huber's psi of a block's outcomes, standardised by the block's MAD
#'
#' A bounded representation: each outcome's departure from the block median,
#' in units of the block's median absolute deviation (MAD, scaled by 1.4826 as
#' [stats::mad()] does, so that it equals the standard deviation for normal
#' outcomes), held within `[-k, k]`. With `k = 1.345`, Huber's M-estimate of
#' location is 95 percent as efficient as the mean when the outcomes are normal.
#'
#' Standardising within block keeps the representation informative whatever
#' the outcome's location and units, and the cap keeps one extreme unit from
#' dominating the block. When more than half of a block's outcomes tie, the MAD
#' is zero; the departures are then divided by the mean absolute departure from
#' the median instead. A block with no spread at all gets zeros.
#'
#' This is the sixth representation of the default set, following Bowers and
#' Burton, "A More Powerful Test for Randomized Experiments" (working paper),
#' and it reproduces that paper's computation.
#'
#' @inheritParams riposte_mean_dist
#' @param k the cap, in MAD units. Huber's 1.345 by default.
#' @return numeric vector the same length as `y`.
#' @export
riposte_huber <- function(y, k = 1.345) {
  ## an unused block level reaches every representation as numeric(0), and
  ## median(numeric(0)) is NA, which the comparisons below cannot handle
  if (length(y) < 1L) return(numeric(0))
  center <- stats::median(y)
  scale <- stats::mad(y, center = center)
  ## the MAD is zero when more than half the block ties at the median; the
  ## mean absolute departure is still positive unless every outcome ties
  if (scale == 0) scale <- mean(abs(y - center))
  if (scale == 0) return(numeric(length(y)))
  pmin(pmax((y - center) / scale, -k), k)
}

#' The default set of outcome representations
#'
#' Returns the six representations of Bowers and Burton, "A More Powerful Test
#' for Randomized Experiments" (working paper), as a named list of functions,
#' each mapping a block's outcome vector to a per-unit score vector of the same
#' length:
#' \describe{
#'   \item{raw}{the outcome itself (identity)}
#'   \item{rank}{within-block mid-ranks, [riposte_rank()]}
#'   \item{mean_dist}{mean L1 distance, [riposte_mean_dist()]}
#'   \item{mean_rank_dist}{mean L1 distance on ranks, [riposte_mean_rank_dist()]}
#'   \item{max_dist}{max distance to the block extremes, [riposte_max_dist()]}
#'   \item{huber}{Huber's psi of the MAD-standardised outcome, [riposte_huber()]}
#' }
#'
#' Pass your own list (or this one extended) to [riposte_test()] or
#' [riposte_components()] to combine a different set. Each element must be a
#' function of one numeric vector returning a numeric vector of the same length;
#' it is applied within each block.
#'
#' @return a named list of representation functions.
#' @seealso [riposte_reps_proposal()] for the earlier set with raw `tanh(y)`.
#' @export
riposte_reps_default <- function() {
  list(
    raw            = function(y) y,
    rank           = riposte_rank,
    mean_dist      = riposte_mean_dist,
    mean_rank_dist = riposte_mean_rank_dist,
    max_dist       = riposte_max_dist,
    huber          = riposte_huber
  )
}

#' The proposal's representation set (raw tanh)
#'
#' [riposte_reps_default()] with Huber's psi replaced by the raw `tanh(y)` that
#' the "Power for Policy" proposal used (as do `manytestsr`'s `pIndepDist()` and
#' `pCombCauchyDist()`), for reproducing the proposal's tables exactly. Prefer
#' [riposte_reps_default()] for real outcomes, whose scale `tanh(y)` does not
#' respect.
#'
#' @return a named list of representation functions.
#' @export
riposte_reps_proposal <- function() {
  reps <- riposte_reps_default()
  ## swap, not append: the proposal's set has six representations too
  reps$huber <- NULL
  reps$tanh <- function(y) tanh(y)
  reps
}

#' Validate a list of outcome representations
#'
#' Checks that `representations` is a non-empty, uniquely named list of
#' functions, each taking one argument. Errors otherwise. The per-block
#' output-length check happens when the functions are applied, in
#' [riposte_score_matrix()].
#'
#' @param representations a named list of functions.
#' @return `representations`, invisibly, if valid.
#' @keywords internal
#' @noRd
riposte_validate_reps <- function(representations) {
  if (!is.list(representations) || length(representations) == 0L)
    stop("`representations` must be a non-empty named list of functions.", call. = FALSE)
  nm <- names(representations)
  if (is.null(nm) || any(nm == "") || anyDuplicated(nm))
    stop("`representations` must have unique, non-empty names.", call. = FALSE)
  if (!all(vapply(representations, is.function, logical(1))))
    stop("every element of `representations` must be a function.", call. = FALSE)
  invisible(representations)
}
