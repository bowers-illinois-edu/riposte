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

#' Two-sided mid-p permutation p-values
#'
#' Given a statistic evaluated across the observed assignment and its
#' re-randomizations (a single vector `v`, length 1 + number of draws), returns
#' the two-sided mid-p permutation p-value of each entry relative to the whole
#' set. The most extreme value (largest `|v|`) gets `0.5 / n`; the least extreme
#' gets `(n - 0.5) / n`. Every value is strictly inside (0, 1), so the Cauchy
#' transform never reaches a pole.
#'
#' The mid-p value is the ordinary upper-tail count minus half the mass at the
#' observed value: `(# strictly more extreme) + 0.5 * (# tied, including self)`,
#' divided by `n`. With `ties.method = "average"` the rank-based form below gives
#' exactly that. The linear statistics riposte combines are continuous, so ties
#' are negligible in practice.
#'
#' @param v numeric vector of a statistic across assignments (observed plus
#'   draws). Two-sidedness is by `abs(v)`.
#' @return numeric vector of mid-p values, the same length as `v`, all in (0, 1).
#' @export
riposte_midp <- function(v) {
  a <- abs(v)
  n <- length(a)
  ## rank by |v| ascending; the largest |v| has the highest rank, hence smallest
  ## p. ties.method = "average" splits tied mass evenly, giving the mid-p.
  (n - rank(a, ties.method = "average") + 0.5) / n
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
#' # The final Cauchy tail is an approximation, not a permutation p-value.
#' @export
riposte_acat_term <- function(p) {
  out <- tan((0.5 - p) * pi)
  small <- p < 1e-4
  out[small] <- 1 / (p[small] * pi)
  big <- p > 1 - 1e-4
  out[big] <- -1 / ((1 - p[big]) * pi)
  out
}
