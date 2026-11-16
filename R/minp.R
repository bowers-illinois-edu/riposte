## minp.R --- the min-p combination and its settings.
##
## THE POINT. The min-p takes the smallest of K p-values at every assignment,
## the observed one and each redrawn one, and its p-value is the share of
## assignments whose smallest p-value is at most the observed smallest. The
## smallest is chosen again inside every assignment, so the min-p rejects a true
## sharp null hypothesis at level alpha no more than a share alpha of the time,
## whatever the K p-values are and however they are computed. The paper (Bowers
## and Burton) recommends the min-p over seven p-values: each of the six
## representations alone and the quadratic combination of the six.
##
## Two settings change which test the min-p is, and riposte_minp_control() holds
## both (Jake's decisions of 2026-11-16):
##   pvalue     how each of the K p-values is computed at every assignment.
##              "midp": its mid-p value among the assignments (riposte's choice
##              since 0.0.0.9009). "randomization": its ordinary randomization
##              p-value, the share of assignments at which it is at least as
##              extreme, ties counted. "asymptotic": the normal tail of the
##              standardized statistic, and the chi-square tail of the quadratic
##              combination, as the paper's code computes them. The level is
##              exact under all three. Only "randomization" makes the min-p's
##              p-value at most K times the smallest of the K, so that the
##              min-p rejects at alpha whenever one of them is at or below
##              alpha / K; under "midp" the bound is 2K times the smallest mid-p
##              value (tests/testthat/test-minp-rules.R).
##   quadratic  whether the quadratic combination is one of the K. It squares
##              each representation's departure, so it responds to departures in
##              both directions and has no place in a one-sided min-p.
##
## Every rule measures a treated sum's departure from the null mean it is given:
## zero without covariance adjustment, and the pooled mean over the observed and
## redrawn assignments with it, as the max and the quadratic do.

#' Settings for the min-p combination
#'
#' The min-p combination (`riposte_test(statistic = "minp")`, the default) takes
#' the smallest of several p-values at the observed assignment and at every
#' re-randomization, and reports the share of assignments whose smallest p-value
#' is at most the observed one. This function sets which p-values it combines
#' and how each is computed.
#'
#' @param pvalue how each p-value inside the min-p is computed at every
#'   assignment. `"midp"` (default): the mid-p value of the statistic among the
#'   assignments, the share more extreme plus half the share tied.
#'   `"randomization"`: the ordinary randomization p-value, the share of
#'   assignments at which the statistic is at least as extreme, ties counted in
#'   full. `"asymptotic"`: the two-sided normal tail of each standardized
#'   statistic (one-sided under a one-sided `alternative`) and the chi-square
#'   tail of the quadratic combination, as in the paper's code. Under every rule
#'   the min-p's own p-value comes from the re-randomizations, so the test keeps
#'   its level. With `"randomization"` the min-p's p-value is also at most K
#'   times the smallest of its K p-values, so the min-p rejects at level alpha
#'   whenever one of them is at or below alpha / K. With `"midp"` that bound
#'   becomes 2K times the smallest mid-p value.
#' @param quadratic whether the quadratic combination's p-value is one of the
#'   p-values combined. `NULL` (default) includes it for a two-sided test and
#'   leaves it out for a one-sided test. The quadratic combination squares each
#'   representation's departure from its re-randomization mean, so it responds
#'   to departures in both directions; asking for it with a one-sided
#'   `alternative` is an error.
#' @return a list of class `riposte_minp_control` with elements `pvalue` and
#'   `quadratic`, for the `minp` argument of [riposte_test()].
#' @examples
#' riposte_minp_control()
#' riposte_minp_control(pvalue = "randomization")
#' riposte_minp_control(quadratic = FALSE)
#' @export
riposte_minp_control <- function(pvalue = c("midp", "randomization", "asymptotic"),
                                 quadratic = NULL) {
  pvalue <- match.arg(pvalue)
  if (!is.null(quadratic) &&
      !(is.logical(quadratic) && length(quadratic) == 1L && !is.na(quadratic)))
    stop("`quadratic` must be NULL, TRUE, or FALSE.", call. = FALSE)
  structure(list(pvalue = pvalue, quadratic = quadratic),
            class = "riposte_minp_control")
}

## The error for a one-sided min-p that asks for the quadratic combination. The
## message says why and names both ways out (Jake, 2026-11-16).
riposte_minp_one_sided_error <- function(alternative) {
  moved <- if (alternative == "greater") "fell" else "rose"
  stop("The min-p cannot include the quadratic combination when alternative = \"",
       alternative, "\". The quadratic combination squares each representation's ",
       "departure from its re-randomization mean, so it responds to departures in ",
       "both directions. A min-p that included it could reject because the treated ",
       "sums ", moved, ". A test of \"", alternative, "\" should not reject for that ",
       "reason. Use quadratic = FALSE (the default for a one-sided test) for a ",
       "one-sided min-p, or alternative = \"two.sided\" to include the quadratic ",
       "combination.", call. = FALSE)
}

## ---- the combination on a statistics matrix ----------------------------------

## Tmat is (1 + nresample) x representations, row 1 the observed assignment. mu
## and Sigma are the null mean and covariance of the treated sums: closed-form
## without adjustment, pooled over the rows with it. The K p-values are kept on
## the log scale, so that a formula p-value too small to represent still orders
## the assignments, and the smallest is compared across rows with riposte's
## relative tie tolerance.
riposte_minp_from_T <- function(Tmat, alternative = "two.sided", pvalue = "midp",
                                quadratic = FALSE, mu = NULL, Sigma = NULL) {
  if (quadratic && alternative != "two.sided")
    riposte_minp_one_sided_error(alternative)
  if (is.null(mu)) mu <- numeric(ncol(Tmat))
  if (is.null(Sigma) && (quadratic || pvalue == "asymptotic"))
    stop("the quadratic combination and the asymptotic rule need the null ",
         "covariance of the treated sums.", call. = FALSE)
  D <- sweep(Tmat, 2, mu)
  logP <- switch(
    pvalue,
    midp = log(apply(D, 2, riposte_midp, alternative = alternative)),
    randomization = log(apply(D, 2, riposte_randp, alternative = alternative)),
    asymptotic = {
      sdv <- sqrt(pmax(diag(Sigma), 0))
      sdv[sdv == 0] <- Inf                 # a sum with no variance never departs
      z <- sweep(D, 2, sdv, "/")
      if (alternative == "two.sided")
        stats::pchisq(z^2, df = 1, lower.tail = FALSE, log.p = TRUE)
      else
        stats::pnorm(if (alternative == "greater") z else -z,
                     lower.tail = FALSE, log.p = TRUE)
    })
  logP <- matrix(logP, nrow = nrow(Tmat), dimnames = list(NULL, colnames(Tmat)))
  df <- NA_integer_
  if (quadratic) {
    ## the same correlation-scale inverse as riposte_quadratic_from_T()
    inv <- riposte_std_pinv(Sigma)
    cen <- sweep(D, 2, inv$sd, "/")
    q <- rowSums((cen %*% inv$R_pinv) * cen)
    df <- inv$rank
    logq <- switch(
      pvalue,
      midp = log(riposte_midp(q, alternative = "greater")),
      randomization = log(riposte_randp(q, alternative = "greater")),
      asymptotic = stats::pchisq(q, df = df, lower.tail = FALSE, log.p = TRUE))
    logP <- cbind(logP, quadratic = logq)
  }
  smallest <- apply(logP, 1, min)
  ## riposte_perm_pvalue() counts values at least as large as the observed one,
  ## so negate: a smaller smallest p-value is more extreme
  list(statistic = exp(smallest[[1]]), p.value = riposte_perm_pvalue(-smallest),
       which = colnames(logP)[which.min(logP[1, ])],
       component_p = stats::setNames(exp(logP[1, ]), colnames(logP)),
       pvalue = pvalue, quadratic = quadratic, df = df)
}

#' Single-step min-p combination
#'
#' @inheritParams riposte_quadratic
#' @param alternative `"two.sided"`, `"greater"`, or `"less"`: the direction of
#'   each representation's p-value.
#' @param pvalue,quadratic as in [riposte_minp_control()], except that
#'   `quadratic` must be `TRUE` or `FALSE` here.
#' @return a list with the observed `statistic` (the smallest p-value), its
#'   permutation `p.value`, `which` (the p-value that was smallest),
#'   `component_p` (each observed p-value), `pvalue`, `quadratic`, and `df` (the
#'   quadratic's rank, `NA` without it).
#' @keywords internal
#' @noRd
riposte_minp <- function(scores, z, block, moments = NULL, draws = NULL,
                         nresample = 1999L, alternative = "two.sided",
                         pvalue = "midp", quadratic = FALSE) {
  scores <- as.matrix(scores)
  if (ncol(scores) == 0L)
    return(list(statistic = NA_real_, p.value = NA_real_, which = NA_character_,
                component_p = numeric(0)))
  moments <- moments %||% riposte_sw_moments(scores, z, block)
  Tmat <- riposte_perm_stats(scores, z, block, nresample, draws)$stats
  riposte_minp_from_T(Tmat, alternative, pvalue, quadratic, moments$mu, moments$Sigma)
}
