## Large-sample calculations for fixed, unadjusted outcome representations.
## The assignment design still determines the exact mean and covariance of
## the treated-score sums. Only the reference used for the p-value changes.

riposte_test_asymptotic <- function(des, sm, statistic) {
  moments <- riposte_sw_moments(sm$scores, des$z, des$block)
  riposte_assert_testable(moments$Sigma)

  ## A score can vary only in a block with no randomized treatment. It then
  ## has zero randomization variance even though its score column is nonzero.
  ## Such scores cannot be standardized and add no evidence to either test.
  keep <- diag(moments$Sigma) > 0
  dropped <- c(sm$dropped, sm$kept[!keep])
  scores <- sm$scores[, keep, drop = FALSE]
  moments$mu <- moments$mu[keep]
  moments$Sigma <- moments$Sigma[keep, keep, drop = FALSE]
  centered <- as.numeric(crossprod(des$z, scores)) - moments$mu

  if (statistic == "quadratic") {
    ## Use the same numerical cutoff as MASS::ginv(), the permutation
    ## quadratic's inverse. The retained directions determine BOTH the
    ## statistic and its degrees of freedom; duplicated scores add neither.
    eig <- eigen(moments$Sigma, symmetric = TRUE)
    retained <- eig$values > max(eig$values) * sqrt(.Machine$double.eps)
    coordinates <- as.numeric(crossprod(eig$vectors[, retained, drop = FALSE],
                                       centered))
    q <- sum(coordinates^2 / eig$values[retained])
    df <- sum(retained)
    result <- list(statistic = q,
                   p.value = stats::pchisq(q, df = df, lower.tail = FALSE),
                   df = df)
  } else {
    standardized <- centered / sqrt(diag(moments$Sigma))
    squared <- standardized^2
    ## Calculate both tails directly. Subtraction from 1 would erase a
    ## small tail, turning a finite Cauchy term into an artificial infinity.
    log_p <- stats::pchisq(squared, df = 1, lower.tail = FALSE, log.p = TRUE)
    log_1mp <- stats::pchisq(squared, df = 1, lower.tail = TRUE, log.p = TRUE)
    result <- riposte_cauchy_log_tails(log_p, log_1mp)
    result$df <- NA_integer_
    result$component_p <- stats::setNames(exp(log_p), colnames(scores))
  }

  return(c(result, list(combination = statistic,
                        condition = riposte_screen(moments$Sigma)$condition,
                        screen = NULL, kept = colnames(scores),
                        dropped = dropped)))
}

## Equal-weight Cauchy combination using both log tails of each p-value.
## Keeping the signed tangent terms on a common scale allows a representable
## combined p-value even when the corresponding statistic is too large to
## store. No permutation fallback or artificial p-value clamp is used.
riposte_cauchy_log_tails <- function(log_p, log_1mp) {
  positive <- log_p < log(0.5)
  negative <- log_p > log(0.5)
  signs <- as.numeric(positive) - as.numeric(negative)
  log_tail <- pmin(log_p, log_1mp)
  log_terms <- rep(-Inf, length(log_p))
  nonzero <- signs != 0

  ## cot(pi * p) is positive below 1/2; above 1/2 its magnitude is
  ## cot(pi * (1-p)). Below 1e-8, cot(x) = 1/x to double precision.
  small <- nonzero & log_tail < log(1e-8)
  log_terms[small] <- -log_tail[small] - log(pi)
  regular <- nonzero & !small
  log_terms[regular] <- log(1 / tan(pi * exp(log_tail[regular])))

  ## At an exact p-value of 1 the defining tangent has limit -Inf, and
  ## the combined upper-tail probability is 1. Opposite infinite terms
  ## have no defined sum, so report that case rather than invent a value.
  infinite <- is.infinite(log_terms) & log_terms > 0
  if (any(infinite)) {
    directions <- unique(signs[infinite])
    if (length(directions) > 1L)
      stop("the asymptotic Cauchy combination has both positive and negative ",
           "infinite terms; use engine = \"permute\".", call. = FALSE)
    return(list(statistic = directions * Inf,
                p.value = if (directions > 0) 0 else 1))
  }

  largest <- max(log_terms)
  if (largest == -Inf) return(list(statistic = 0, p.value = 0.5))
  scaled_mean <- mean(signs * exp(log_terms - largest))
  if (scaled_mean == 0) return(list(statistic = 0, p.value = 0.5))
  direction <- sign(scaled_mean)
  log_magnitude <- largest + log(abs(scaled_mean))
  statistic <- direction * exp(log_magnitude)

  ## For large |T|, its small Cauchy tail is 1/(pi*|T|) to double
  ## precision. Calculating that tail from log(|T|) also handles overflow
  ## in the statistic without rounding a representable p-value to zero.
  if (log_magnitude > log(1e8)) {
    log_small_tail <- -log_magnitude - log(pi)
    p <- if (direction > 0) exp(log_small_tail) else -expm1(log_small_tail)
  } else {
    p <- stats::pcauchy(statistic, lower.tail = FALSE)
  }
  return(list(statistic = statistic, p.value = p))
}
