## pvalues.R --- permutation p-value helpers for the combinations.
##
## These helpers support the permutation Cauchy combination. The asymptotic
## engine uses the separate log-tail calculation in asymptotic-engine.R:
##   - riposte_midp(): the two-sided mid-p permutation p-value, strictly inside
##     (0, 1). The Cauchy transform tan((0.5 - p) * pi) has poles at 0 and 1; a
##     plain permutation p-value can hit exactly 1 (the least extreme component),
##     which a naive clamp turns into a huge negative term that swamps the
##     average. The mid-p value never reaches a pole, so no clamp is needed. It
##     holds the nominal level and matches the clamp's power (verified in the
##     proposal sims).
##   - riposte_acat_term(): the pole-aware Cauchy transform, for the rare case of
##     combining p-values that did reach 0 or 1 (e.g. p-values from outside
##     riposte). With mid-p inputs it almost never triggers; kept for robustness.

#' Mid-p permutation p-values
#'
#' Given a statistic evaluated across the observed assignment and its
#' re-randomizations (a single vector `v`, length 1 + number of draws), returns
#' the mid-p permutation p-value of each entry relative to the whole set:
#' two-sided by default, or one-sided with `alternative = "greater"` (large
#' values are extreme) or `"less"` (small values are extreme). The most extreme value (largest `|v|`) gets `0.5 / n`; the least extreme
#' gets `(n - 0.5) / n`. Every value is strictly inside (0, 1), so the Cauchy
#' transform never reaches a pole.
#'
#' The mid-p value is the ordinary upper-tail count minus half the mass at the
#' observed value: `(# strictly more extreme) + 0.5 * (# tied, including self)`,
#' divided by `n`. Values within about `1.5e-8` of each other, relative to
#' their size when it exceeds one, count as tied, because rank scores make many
#' re-randomizations tie exactly, and floating-point summation in a different
#' order can separate such ties in their last bits.
#'
#' @param v numeric vector of a statistic across assignments (observed plus
#'   draws). Two-sidedness is by `abs(v)`.
#' @param alternative `"two.sided"` (default), `"greater"`, or `"less"`.
#' @return numeric vector of mid-p values, the same length as `v`, all in (0, 1).
#'   Without ties, the `"greater"` and `"less"` values add to one.
#' @export
riposte_midp <- function(v, alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  a <- switch(alternative, two.sided = abs(v), greater = v, less = -v)
  n <- length(a)
  ## mid-p: (# strictly more extreme + half the # tied, self included) / n, with
  ## values within floating-point tolerance of a_i treated as tied with it
  tol <- riposte_tie_tol(a)
  sorted <- sort(a)
  more <- n - findInterval(a + tol, sorted)                      # a_j > a_i + tol
  at_least <- n - findInterval(a - tol, sorted, left.open = TRUE) # a_j >= a_i - tol
  (more + 0.5 * (at_least - more)) / n
}

## Tolerance within which another statistic counts as tied with x. Rank scores
## make many re-randomizations tie exactly in exact arithmetic (an assignment and
## its within-block complement give treated sums T and -T, for instance), and
## floating-point arithmetic separates such ties in their last bits; without a
## tolerance that noise decides whether a tied re-randomization counts as at
## least as extreme. The tolerance is sqrt(machine epsilon), about 1.5e-8,
## relative to |x| (absolute below 1), as all.equal() uses: far above rounding
## noise, far below any difference moving a unit between groups can make. The
## tolerance belongs at the comparison rather than in the sums, because the
## quadratic form, the standardized max, and the Cauchy average all round after
## the sums are formed.
riposte_tie_tol <- function(x) {
  sqrt(.Machine$double.eps) * pmax(1, abs(x))
}

#' Pole-aware Cauchy (ACAT) transform of p-values
#'
#' `tan((0.5 - p) * pi)`, with the accurate small-/large-p limits substituted
#' near the poles (this is what Liu et al.'s ACAT implementation does):
#' `1 / (p * pi)` as `p -> 0` and `-1 / ((1 - p) * pi)` as `p -> 1`. With mid-p
#' inputs the substitution almost never fires; it makes the transform safe if
#' riposte is ever asked to combine p-values produced elsewhere that reach 0 or 1.
#'
#' @param p numeric vector of p-values.
#' @return numeric vector of Cauchy-transformed values, the same length as `p`.
#' @examples
#' # The paper's Hybrid test gives equal weight to seven p-values:
#' # six individual representations and their quadratic combination.
#' dat <- data.frame(
#'   outcome = c(-3, -1, 0, 0.5, 2, 4, 7, 9,
#'               -2, -1, 0, 1, 2, 3, 4, 5, 8, 13),
#'   treated = c(1, 0, 1, 0, 0, 1, 0, 0,
#'               1, 0, 1, 1, 0, 1, 0, 1, 1, 0),
#'   block = factor(rep(1:2, c(8, 10)))
#' )
#' six <- riposte_test(outcome ~ treated | block, dat,
#'                     statistic = "cauchy", engine = "asymptotic")
#' quad <- riposte_test(outcome ~ treated | block, dat,
#'                      statistic = "quadratic", engine = "asymptotic")
#' hybrid_inputs <- c(six$component_p, quadratic = quad$p.value)
#' hybrid_statistic <- mean(riposte_acat_term(hybrid_inputs))
#' pcauchy(hybrid_statistic, lower.tail = FALSE)
#'
#' # These are seven inputs, not the two combined p-values with equal weight.
#' # The final Cauchy tail is an approximation, not a permutation p-value, and a
#' # single input of exactly 1 makes it 1. riposte_test(statistic = "hybrid")
#' # computes the Hybrid either from re-randomization or, with
#' # engine = "asymptotic", with the truncated conversion of
#' # riposte_truncated_cauchy(), which keeps a p-value of 1 finite.
#' @export
riposte_acat_term <- function(p) {
  out <- tan((0.5 - p) * pi)
  small <- p < 1e-4
  out[small] <- 1 / (p[small] * pi)
  big <- p > 1 - 1e-4
  out[big] <- -1 / ((1 - p[big]) * pi)
  out
}
