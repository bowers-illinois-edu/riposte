## moments.R --- closed-form Strasser-Weber permutation moments.   [STUB]
##
## THE POINT. The quadratic omnibus and the screen need the permutation mean and
## covariance of the representation linear statistics. The package must compute
## these in CLOSED FORM (Strasser & Weber 1999), as functions of the centred
## scores and the design (block sizes and treated counts), NOT by Monte Carlo.
## Closed-form moments make the quadratic exactly exact --- the proposal's
## Monte-Carlo moments make it anti-conservative by O(1/nresample) --- and make
## the screen's condition number exactly ancillary to the realized assignment.
##
## For a single block with centred scores s (column-centred so sum_i s_i = 0),
## treated count m, block size k, the permutation mean of sum_{treated} s is 0
## and the permutation covariance between representations a, b is
##   Cov = m (k - m) / (k - 1) * (1/k) * sum_i s_{a,i} s_{b,i}
## (the finite-population covariance of a sum over a simple random sample without
## replacement of size m from k items with values s, summed across the joint).
## Blocks are independent under the design, so the design-level mu and Sigma are
## the sums of the per-block contributions. Implement and verify against coin's
## expectation()/covariance() in the tests.
##
## DECISION PENDING (see the checkpoint notes): implement these formulas natively
## (coin -> Suggests, used only for the equivalence test) versus extract them from
## a coin IndependenceTest at runtime (coin -> Imports). The native route keeps
## riposte light and owns its exactness; it is the current plan.

#' Closed-form Strasser-Weber permutation mean and covariance
#'
#' Computes the exact permutation mean vector and covariance matrix of the
#' representation linear statistics `T_c = sum_i z_i * s_{c,i}`, from the scores
#' and the design, without re-randomizing. The permutation fixes each block's
#' treated count, so within block `z` is a simple random sample of size `m_b`
#' from `n_b` units. For that sampling, with `p_b = m_b / n_b`,
#' \itemize{
#'   \item `E[T_c] = sum_b p_b * sum_{i in b} s_{c,i}`, and
#'   \item `Cov(T_a, T_b) = sum_b w_b * sum_{i in b} (s_{a,i} - sbar_a)(s_{b,i} - sbar_b)`,
#'         where `w_b = m_b (n_b - m_b) / (n_b (n_b - 1))` and `sbar` is the
#'         within-block mean.
#' }
#' The covariance uses within-block-centred cross products, so it is the same
#' whether `scores` are already centred (as from [riposte_score_matrix()]) or
#' not; the function centres internally. When `scores` are centred, the mean is
#' zero, as riposte's combinations assume.
#'
#' Verified against `coin`'s `expectation()` and `covariance()` (the same
#' linear-statistic framework) to machine precision; see the tests.
#'
#' @param scores numeric matrix of scores (rows = units, columns =
#'   representations), e.g. `riposte_score_matrix()$scores`.
#' @param z 0/1 treatment vector (its per-block sums give the treated counts).
#' @param block factor of block labels.
#' @return a list with `mu` (length = number of representations) and `Sigma`
#'   (representation-by-representation covariance matrix).
#' @export
riposte_sw_moments <- function(scores, z, block) {
  scores <- as.matrix(scores)
  storage.mode(scores) <- "double"
  z <- as.numeric(z)
  block <- as.factor(block)
  n <- nrow(scores)
  if (length(z) != n || length(block) != n)
    stop("`scores`, `z`, and `block` must describe the same units.", call. = FALSE)

  cols <- colnames(scores)
  C <- ncol(scores)
  mu <- numeric(C)
  Sigma <- matrix(0, C, C)

  ## blocks are independent under the design, so mu and Sigma are sums of
  ## per-block contributions
  for (ix in split(seq_len(n), block)) {
    nb <- length(ix)
    mb <- sum(z[ix])
    Sb <- scores[ix, , drop = FALSE]
    bmean <- colMeans(Sb)
    mu <- mu + (mb / nb) * colSums(Sb)         # conditional mean within block
    if (nb < 2L) next                          # a singleton has no variation
    w <- mb * (nb - mb) / (nb * (nb - 1))       # 0 if all treated or all control
    if (w == 0) next
    Sc <- sweep(Sb, 2, bmean)                   # centre within block
    Sigma <- Sigma + w * crossprod(Sc)
  }

  names(mu) <- cols
  dimnames(Sigma) <- list(cols, cols)
  list(mu = mu, Sigma = Sigma)
}
