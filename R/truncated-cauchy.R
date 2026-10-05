## truncated-cauchy.R --- the truncated Cauchy combination of Gui, Jiang, and
## Wang (2025).
##
## THE POINT. Liu and Xie's Cauchy combination converts each p-value p to
## tan((0.5 - p) * pi) and averages. That conversion makes one small p-value
## dominate the average, which is the reason for using it, but it treats a
## p-value near 1 the same way with the sign reversed: 0.999 converts to about
## -318 and 1 to minus infinity, so one representation whose treated sum sits
## at its re-randomization mean sets the combined p-value to 1 whatever the
## others say. On coarse outcomes with few units such a sum is common, and its
## large-sample p-value is then exactly 1.
##
## Gui, Jiang, and Wang (2025, Biometrika 112(4), asaf038, Section 2.2) convert
## p instead with the quantile function of a Cauchy distribution cut off below
## its (1 - t) quantile, tan((0.5 - t) * pi). With t = 0.9, their
## recommendation, the cut-off point is tan(-0.4 * pi) = -3.08, and that is
## what p = 1 converts to. Solving F(x) = 1 - p for the cut-off distribution's
## CDF F gives the conversion x = tan((0.5 - t * p) * pi). The cut-off
## distribution keeps a share t of the Cauchy's probability, so its upper tail
## is the Cauchy's divided by t, and for n p-values whose converted values sum
## to S the combined p-value is min(1, n * P(Cauchy > S) / t), the heavy-tail
## approximation to the tail of a sum of n such variables (their Section 2.3).
## With one p-value this returns that p-value exactly.

## The truncated combination from the logarithms of the p-values, so that a
## p-value too small to convert directly (its converted value overflows) still
## gives a representable combined p-value. Returns the mean converted value as
## the statistic, as the untruncated Cauchy combination reports its mean.
riposte_truncated_cauchy_log <- function(log_p, truncation = 0.9) {
  n <- length(log_p)
  if (any(log_p == -Inf)) return(list(statistic = Inf, p.value = 0))
  ## u = t * p, at most t < 1, so the conversion cot(pi * u) never reaches the
  ## pole at u = 1; below u = 1e-8, cot(pi * u) = 1 / (pi * u) to double precision
  log_u <- log(truncation) + log_p
  small <- log_u < log(1e-8)
  regular <- 1 / tan(pi * exp(log_u[!small]))
  if (!any(small)) {
    S <- sum(regular)
    p <- n * stats::pcauchy(S, lower.tail = FALSE) / truncation
    return(list(statistic = S / n, p.value = min(1, p)))
  }
  ## at least one converted value exceeds 3e7: add the large ones on the log
  ## scale, and use the Cauchy tail 1 / (pi * S), exact to double precision there
  log_big <- -log_u[small] - log(pi)
  top <- max(log_big)
  log_S <- top + log(sum(exp(log_big - top)) + sum(regular) * exp(-top))
  log_p_comb <- log(n) - log_S - log(pi) - log(truncation)
  return(list(statistic = exp(log_S - log(n)), p.value = min(1, exp(log_p_comb))))
}

#' Truncated Cauchy combination of p-values
#'
#' Combines p-values with the truncated Cauchy conversion of Gui, Jiang, and
#' Wang (2025). Each p-value \eqn{p} is converted to
#' \eqn{\tan((0.5 - t p)\pi)}, where \eqn{t} is `truncation`; with \eqn{n}
#' p-values whose converted values sum to \eqn{S}, the combined p-value is
#' \eqn{\min(1, n P(C > S) / t)}, with \eqn{C} a standard Cauchy variable.
#'
#' Liu and Xie's Cauchy combination converts \eqn{p} to
#' \eqn{\tan((0.5 - p)\pi)}, which is minus infinity at \eqn{p = 1}: one
#' p-value of 1 then makes the combined p-value 1 whatever the others are, and
#' a p-value near 1 does almost as much (0.999 converts to about -318). On
#' outcomes with few distinct values a representation's treated sum can equal
#' its mean exactly, and its large-sample p-value is then 1. With
#' \eqn{t = 0.9}, the value Gui, Jiang, and Wang recommend, \eqn{p = 1}
#' converts to \eqn{\tan(-0.4\pi) = -3.08}. `truncation = 1` gives Liu and
#' Xie's untruncated combination, referred to the standard Cauchy as
#' [riposte_test()] with `engine = "asymptotic"` did before riposte 0.0.0.9008.
#'
#' @param p numeric vector of p-values in \[0, 1\].
#' @param truncation the share \eqn{t} of the Cauchy distribution kept, in
#'   (0, 1\]; 0.9 by default.
#' @return the combined p-value.
#' @references Gui, L., Jiang, Y., and Wang, J. (2025). Aggregating dependent
#'   signals with heavy-tailed combination tests. *Biometrika*, 112(4),
#'   asaf038. \doi{10.1093/biomet/asaf038}
#'
#'   Fang, Y., Chang, C., Park, Y., and Tseng, G. C. (2023). Heavy-tailed
#'   distribution for combining dependent p-values with asymptotic robustness.
#'   *Statistica Sinica*, 33, 1115-1142. \doi{10.5705/ss.202022.0046}
#' @examples
#' # four small p-values and one representation that finds no difference
#' p <- c(0.01, 0.01, 0.01, 0.01, 1)
#' riposte_truncated_cauchy(p)                  # 0.013
#' riposte_truncated_cauchy(p, truncation = 1)  # Liu and Xie's version: 1
#' @export
riposte_truncated_cauchy <- function(p, truncation = 0.9) {
  riposte_check_truncation(truncation)
  if (!is.numeric(p) || anyNA(p) || any(p < 0 | p > 1))
    stop("p must be numeric p-values in [0, 1] with no NA.", call. = FALSE)
  if (truncation == 1)
    return(riposte_cauchy_log_tails(log(p), log1p(-p))$p.value)
  return(riposte_truncated_cauchy_log(log(p), truncation)$p.value)
}

riposte_check_truncation <- function(truncation) {
  if (!is.numeric(truncation) || length(truncation) != 1L || is.na(truncation) ||
      truncation <= 0 || truncation > 1)
    stop("cauchy_truncation (truncation) must be a single number in (0, 1].",
         call. = FALSE)
  invisible(truncation)
}
