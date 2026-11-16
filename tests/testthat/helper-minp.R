## helper-minp.R --- brute-force min-p for the tests of its p-value rules.
##
## The min-p takes the smallest of K p-values at every assignment and refers
## it to the same smallest value at the other assignments. Three rules give
## the K p-values (riposte_minp_control()):
##   midp           the mid-p value of each statistic among the assignments;
##   randomization  the ordinary randomization p-value: the share of the
##                  assignments at which the statistic is at least as extreme,
##                  ties counted, the observed assignment included;
##   asymptotic     the normal tail of each standardized statistic, and the
##                  chi-square tail of the quadratic combination, as the
##                  paper's code computes them.
## These functions compute all of that from a full enumeration with nothing
## from riposte but riposte_midp(), so the tests compare riposte's
## combination with a count. The null mean and covariance come from the
## enumeration itself: over every assignment the design allows, the mean and
## the (divide-by-N) covariance of the treated sums ARE the Strasser-Weber
## moments, so the brute force does not reuse riposte_sw_moments().

## mean and covariance of the treated sums over every assignment
enum_moments <- function(Tmat) {
  mu <- colMeans(Tmat)
  cen <- sweep(Tmat, 2, mu)
  list(mu = mu, Sigma = crossprod(cen) / nrow(Tmat))
}

## share of entries at least as large as each entry, ties counted, self
## included: the ordinary randomization p-value of a statistic for which
## large values are extreme. Ties within 1e-9 (relative) count as ties; riposte
## uses sqrt(.Machine$double.eps), about 1.5e-8. The two agree on the test
## designs, where exact ties (an assignment and its within-block complement)
## differ only by rounding and distinct values differ by far more than 1e-8.
share_at_least <- function(a) {
  tol <- 1e-9 * pmax(1, abs(a))
  vapply(seq_along(a), function(i) mean(a >= a[i] - tol[i]), numeric(1))
}

## The K p-values at every assignment (one row per assignment). With
## quadratic = TRUE a last column holds the quadratic combination's p-value;
## the quadratic form inverts Sigma with solve(), so the test designs keep
## Sigma of full rank.
brute_member_p <- function(Tmat, rule, alternative = "two.sided", quadratic = FALSE) {
  mom <- enum_moments(Tmat)
  z <- sweep(sweep(Tmat, 2, mom$mu), 2, sqrt(diag(mom$Sigma)), "/")
  directed <- switch(alternative, two.sided = abs(z), greater = z, less = -z)
  P <- switch(rule,
    midp = apply(Tmat, 2, riposte_midp, alternative = alternative),
    randomization = apply(directed, 2, share_at_least),
    asymptotic = if (alternative == "two.sided")
      stats::pchisq(z^2, df = 1, lower.tail = FALSE)
    else stats::pnorm(directed, lower.tail = FALSE))
  P <- matrix(P, nrow = nrow(Tmat), dimnames = list(NULL, colnames(Tmat)))
  if (quadratic) {
    stopifnot(qr(mom$Sigma)$rank == ncol(Tmat))
    cen <- sweep(Tmat, 2, mom$mu)
    q <- rowSums((cen %*% solve(mom$Sigma)) * cen)
    pq <- switch(rule,
      midp = riposte_midp(q, alternative = "greater"),
      randomization = share_at_least(q),
      asymptotic = stats::pchisq(q, df = ncol(Tmat), lower.tail = FALSE))
    P <- cbind(P, quadratic = pq)
  }
  P
}

## the min-p's p-value with each assignment taken in turn as the observed one:
## the share of assignments whose smallest p-value is at most that one's
brute_minp_all <- function(P) {
  s <- apply(P, 1, min)
  vapply(s, function(x) mean(s <= x + 1e-9 * max(1, x)), numeric(1))
}

## the draws matrix with assignment j as the observed one (first column)
observed_first <- function(G, j) cbind(G[, j], G[, -j])

## the indicator of being above the block's third quartile: a representation
## whose treated sum takes few values
top_quartile <- function(v) as.numeric(v > stats::quantile(v, 0.75))

## a small block-randomized experiment with an effect, for the riposte_test()
## interface tests. riposte_test(seed = ) switches the generator to
## L'Ecuyer-CMRG, so the design sets the generator itself; otherwise the data
## for a given seed would depend on which tests ran before.
minp_interface_design <- function(seed, B = 4L, nb = 8L, effect = 0.8) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  block <- factor(rep(seq_len(B), each = nb))
  z <- integer(B * nb)
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L
  }
  data.frame(y = rnorm(B * nb) + effect * z, z = z, b = block,
             x = rnorm(B * nb))
}

## ---- the block statistic ---------------------------------------------------

## The block statistic at every assignment (one column of G per assignment),
## written out from its definition: in each block, the treated sum of the
## score's departures from the block mean, squared, divided by its variance
## m_b (n_b - m_b) / (n_b - 1) x sum of squared departures / n_b, and added
## over the blocks whose variance is positive. Returns the statistic and the
## number of blocks in the sum, the chi-square's degrees of freedom.
brute_block_stat <- function(score, block, G) {
  stat <- numeric(ncol(G))
  used <- 0L
  for (b in levels(factor(block))) {
    ix <- which(block == b)
    dev <- score[ix] - mean(score[ix])
    n_b <- length(ix)
    m_b <- sum(G[ix, 1])
    v_b <- m_b * (n_b - m_b) / (n_b - 1) * sum(dev^2) / n_b
    if (v_b <= 0) next
    used <- used + 1L
    for (j in seq_len(ncol(G))) stat[j] <- stat[j] + sum(G[ix, j] * dev)^2 / v_b
  }
  list(stat = stat, df = used)
}

## the block statistic's p-value at every assignment under each rule; large
## values are extreme
brute_block_p <- function(bs, rule) {
  switch(rule,
    midp = riposte_midp(bs$stat, alternative = "greater"),
    randomization = share_at_least(bs$stat),
    asymptotic = stats::pchisq(bs$stat, df = bs$df, lower.tail = FALSE))
}

## within-block ranks, the score the block tests below use
within_block_rank <- function(y, block) stats::ave(y, block, FUN = rank)

## a block-randomized experiment in which treatment raises outcomes in half the
## blocks and lowers them in the rest, so the within-block comparisons cancel
## when added over blocks
help_and_harm_design <- function(seed, B = 10L, nb = 20L, effect = 0.8) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  block <- factor(rep(seq_len(B), each = nb))
  z <- integer(B * nb)
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L
  }
  direction <- ifelse(as.integer(block) <= B %/% 2L, 1, -1)
  data.frame(y = rnorm(B * nb) + effect * direction * z, z = z, b = block)
}
