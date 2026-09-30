## screen.R --- choose quadratic vs Cauchy from the conditioning of Sigma. [STUB]
##
## THE POINT. The quadratic must invert Sigma; the six representations overlap
## (condition number of their correlation matrix ~ 250 in the proposal design),
## so the inverse can be unstable --- exactly when the Cauchy, which inverts
## nothing, is safer. The screen picks between them from kappa(cov2cor(Sigma)).
##
## WHY THIS IS VALID (state it in the docs): the condition number is computed
## from the outcomes (fixed under the sharp null) and the design via the
## CLOSED-FORM Sigma, NOT from the realized assignment, so it is constant across
## the randomization reference --- an ancillary quantity. Choosing the combination
## from an ancillary statistic does not disturb the exact level.
##
##   threshold  hard switch: kappa > threshold -> Cauchy, else quadratic.
##   shrink     smooth: shrink Sigma toward its diagonal (Sigma + lambda I, or
##              Ledoit-Wolf), lambda chosen by the conditioning. lambda = 0 is the
##              quadratic; large lambda behaves like the Cauchy. Preferred default.
##
## OPEN RESEARCH QUESTION (Project 1 deliverable): the threshold / shrinkage rule.
## Ship a sensible default and let the user tune. The default lambda rule chosen
## here is a placeholder for Jake to review, not a settled result.

#' Control parameters for the quadratic-vs-Cauchy screen
#'
#' @param method `"shrink"` (smooth, the default) or `"threshold"` (hard switch).
#' @param threshold condition-number cutoff for `method = "threshold"`.
#' @param lambda optional fixed shrinkage; if `NULL`, chosen from the
#'   conditioning of Sigma.
#' @return a list of screen settings.
#' @export
riposte_screen_control <- function(method = c("shrink", "threshold"),
                                   threshold = 100, lambda = NULL) {
  method <- match.arg(method)
  list(method = method, threshold = threshold, lambda = lambda)
}

#' Choose (or blend) the combination from the conditioning of Sigma
#'
#' Given a permutation covariance `Sigma` of the representations, returns the
#' screen's decision. Under `"threshold"`, if the condition number of the
#' correlation matrix exceeds `threshold`, choose the Cauchy (which inverts
#' nothing), else the quadratic. Under `"shrink"` (the default), shrink the
#' correlation matrix toward the identity, `R_a = (1 - a) R + a I`, choosing the
#' intensity `a` so the condition number drops to `threshold`; `a = 0` is the full
#' quadratic and `a = 1` is the diagonal covariance (which ignores the
#' cross-representation dependence, as the Cauchy does). The shrunk covariance is
#' returned for the quadratic to use.
#'
#' Choosing the combination from the condition number does not disturb the test's
#' exact level, but the reason differs by path. Without covariance adjustment,
#' `Sigma` is the closed-form permutation covariance --- a function of the outcomes
#' (fixed under the sharp null) and the design, NOT of the realized assignment ---
#' so the condition number is constant across the randomization reference, an
#' ancillary quantity. With covariance adjustment, `Sigma` is the pooled empirical
#' covariance of the residual statistics, which depends on the realized draws; it
#' is no longer ancillary, but it is a symmetric function of the exchangeable
#' pooled set, which is the property exactness actually needs (see
#' `dev/theory-adjusted-exactness.md`).
#'
#' The shrinkage-intensity rule here is a sensible default, not the settled
#' Project 1 result; expect it to change. Set `lambda` in
#' [riposte_screen_control()] to fix the intensity yourself. Exact level holds for
#' any shrinkage rule; the rule affects power only.
#'
#' @param Sigma a permutation covariance of the representations (the closed-form
#'   covariance without adjustment, the pooled empirical covariance with it).
#' @param control a [riposte_screen_control()] list.
#' @return a list with `method`, the `condition` number of the correlation
#'   matrix, and either `choice` (`"quadratic"`/`"cauchy"`, for `"threshold"`) or
#'   the shrinkage `lambda`, the shrunk `Sigma`, and the achieved
#'   `condition_shrunk` (for `"shrink"`).
#' @export
riposte_screen <- function(Sigma, control = riposte_screen_control()) {
  Sigma <- as.matrix(Sigma)
  R <- stats::cov2cor(Sigma)
  ev <- pmax(eigen(R, symmetric = TRUE, only.values = TRUE)$values, 0)
  lmax <- max(ev); lmin <- min(ev)
  condition <- if (lmin <= 0) Inf else lmax / lmin

  if (control$method == "threshold") {
    choice <- if (condition > control$threshold) "cauchy" else "quadratic"
    return(list(method = "threshold", condition = condition, choice = choice,
                Sigma = Sigma, lambda = 0))
  }

  ## shrink toward the identity in correlation space until the condition number
  ## reaches the target. Eigenvalues of R_a are (1 - a) ev + a, so the condition
  ## number ((1-a)lmax + a)/((1-a)lmin + a) = target solves linearly for a.
  target <- control$threshold
  if (!is.null(control$lambda)) {
    a <- min(max(control$lambda, 0), 1)
  } else if (!is.finite(condition) || condition > target) {
    denom <- (1 - lmax) - target * (1 - lmin)
    a <- if (denom == 0) 1 else (target * lmin - lmax) / denom
    a <- min(max(a, 0), 1)
  } else {
    a <- 0
  }

  Ra <- (1 - a) * R + a * diag(nrow(R))
  D <- sqrt(diag(Sigma))
  Sigma_shrunk <- outer(D, D) * Ra
  dimnames(Sigma_shrunk) <- dimnames(Sigma)
  ev_a <- pmax(eigen(Ra, symmetric = TRUE, only.values = TRUE)$values, 0)
  condition_shrunk <- if (min(ev_a) <= 0) Inf else max(ev_a) / min(ev_a)

  list(method = "shrink", condition = condition, lambda = a,
       Sigma = Sigma_shrunk, condition_shrunk = condition_shrunk)
}
