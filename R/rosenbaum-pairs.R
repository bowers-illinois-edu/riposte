## rosenbaum-pairs.R
## Rosenbaum (2012) matched-pair tests, implemented in riposte.
##
## Two procedures for matched pairs, both built on signed statistics
## T = sum_i B_i h_i, where B_i = 1{y_i > 0} and h_i is a score that depends on the
## rank of |y_i|:
##   (a) "testing one hypothesis twice" (Rosenbaum 2012, Biometrika 99:763-774):
##       standardize L statistics and refer their MAXIMUM to its distribution when
##       every sign is positive with chance gamma / (1 + gamma), either by listing or
##       drawing sign patterns (riposte_pairs_perm_test) or, as in the paper and the
##       CRAN package `testtwice`, by a multivariate normal orthant probability
##       (riposte_pairs_bound).
##   (b) the exact adaptive test (Rosenbaum 2012, Ann. Appl. Stat. 6:83-105): reject if
##       Noether's B1 >= kB or Brown's T = 2*B1 + B2 >= kT, with (kB, kT) from exact
##       binomial arithmetic. Reference: the paper's Section 4.2 worked example.
## We cite `testtwice` and the papers and re-implement in riposte; `testtwice` is only
## used (optionally) in the tests as an independent check.
##
## Matched pairs only: a pair's treated-minus-control difference y_i, with
## sensitivity parameter gamma >= 1 (gamma = 1 is the randomization test). This is a
## different object from riposte's block-randomized assignment contrasts, so it lives
## in its own set of functions (riposte_pairs_*), not in riposte_score_matrix().
## Zero differences get score 0 (they carry no sign information), as in `testtwice`.

## Shared argument checks for the public matched-pair entry points

#' @keywords internal
#' @noRd
riposte_check_gamma <- function(gamma) {
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 1)
    stop("`gamma` must be a single finite number >= 1.", call. = FALSE)
  invisible(gamma)
}

#' @keywords internal
#' @noRd
riposte_check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single number strictly between 0 and 1.", call. = FALSE)
  invisible(alpha)
}
#' @keywords internal
#' @noRd
riposte_check_lambda <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 2L || anyNA(lambda) ||
      !(lambda[1] > 0 && lambda[1] < lambda[2] && lambda[2] <= 1))
    stop("`lambda` must be two numbers with 0 < lambda[1] < lambda[2] <= 1.", call. = FALSE)
  invisible(lambda)
}
## Is rank rk at or above the fraction f of the I pairs (rk >= f * I)? Brown's and
## Noether's scores and the adaptive test's bands all use this rule, so that T = 2 B1 +
## B2 stays twice the sum of Brown's scores. Ranks are multiples of 1/2, so the small
## tolerance only settles a rank lying exactly on f * I (rank 20 of 30 at f = 2/3), which
## floating point would otherwise decide by how f * I happens to round.
#' @keywords internal
#' @noRd
riposte_rank_at_least <- function(rk, I, f) rk >= f * I - 1e-8

#' Signed-rank-type scores for matched-pair differences
#'
#' Scores `h_i` that depend on the rank of `|y_i|`, for use with
#' [riposte_pairs_bound()] and [riposte_pairs_test()].
#'
#' @param y numeric vector of treated-minus-control matched-pair differences.
#' @param type `"wilcoxon"` (ranks of `|y|`), `"brown"` (Brown / Markowski-
#'   Hettmansperger combined quantile average, quantiles 1/3 and 2/3), `"noether"`
#'   (0/1 for ranks in the top third), or `"ustat"` (Rosenbaum's 2011 U-statistic,
#'   needs `m`, `m1`, `m2`).
#' @param m,m1,m2 U-statistic indices with `1 <= m1 <= m2 <= m`; `(8, 7, 8)` and
#'   `(8, 6, 7)` are the two highlighted in Rosenbaum (2012).
#' @param q1,q2 quantiles for `"brown"` (defaults 1/3, 2/3) and `"noether"` (`q1`
#'   only, default 2/3).
#' @param exact for `"ustat"`, use the integer form `sum choose(a-1, l-1) *
#'   choose(I-a, m-l)` of the paper (`TRUE`) or the smooth approximation (`FALSE`,
#'   what `testtwice::multrnk()` does by default).
#' @return a numeric vector of scores, scaled to have maximum 1 (the standardized
#'   statistics do not depend on the scale), with 0 for zero differences.
#' @export
riposte_pairs_scores <- function(y, type = c("wilcoxon", "brown", "noether", "ustat"),
                                 m = 8, m1 = 7, m2 = 8, q1 = NULL, q2 = NULL,
                                 exact = TRUE) {
  type <- match.arg(type)
  if (!is.numeric(y) || anyNA(y)) stop("`y` must be numeric without missing values.", call. = FALSE)
  I <- length(y)
  rk <- rank(abs(y))
  nz <- as.numeric(abs(y) > 0)
  h <- switch(type,
              wilcoxon = rk,
              brown = {
                a <- q1 %||% (1 / 3); b <- q2 %||% (2 / 3)
                riposte_rank_at_least(rk, I, a) + riposte_rank_at_least(rk, I, b)
              },
              noether = {
                a <- q1 %||% (2 / 3); riposte_rank_at_least(rk, I, a) * 1
              },
              ustat = {
                if (!(m1 >= 1 && m1 <= m2 && m2 <= m && all(c(m, m1, m2) == round(c(m, m1, m2)))))
                  stop("need integers 1 <= m1 <= m2 <= m.", call. = FALSE)
                if (exact) {
                  out <- numeric(I)
                  for (l in m1:m2)
                    out <- out + exp(lchoose(rk - 1, l - 1) + lchoose(I - rk, m - l))
                  out
                } else {
                  p <- rk / I; out <- numeric(I)
                  for (l in m1:m2)
                    out <- out + l * choose(m, l) * p^(l - 1) * (1 - p)^(m - l)
                  out
                }
              })
  h <- h * nz
  if (max(h) > 0) h <- h / max(h)
  h
}

#' Maximum-of-L sensitivity bound (Rosenbaum 2012, procedure 6)
#'
#' Given L score columns `H` for the same matched pairs, computes each statistic
#' `T_l = sum B_i H[i, l]`, its standardized worst-case deviate at sensitivity
#' `gamma`, and the upper bound on the one-sided p-value of their maximum from the
#' multivariate normal orthant probability implied by `crossprod(H)`.
#'
#' This is a normal approximation (as in the paper and in `testtwice::tt()`).
#' [riposte_pairs_perm_test()] computes the same bound from the distribution of the
#' signs instead, which is what [riposte_test()] does unless `engine = "asymptotic"`.
#'
#' @param y matched-pair differences.
#' @param H numeric matrix, one row per pair, one column per statistic (`2 <= L`).
#' @param gamma sensitivity parameter, `>= 1`.
#' @param algorithm an `mvtnorm` algorithm object. Default: `mvtnorm::Miwa()` for
#'   `L <= 5` (what `testtwice` uses), else `mvtnorm::GenzBretz()` (Miwa is
#'   exponentially slow in L; GenzBretz is Monte Carlo, so set a seed).
#' @return a list with `p.value`, `deviates`, and `cor`.
#' @references Rosenbaum, P. R. (2012) Testing one hypothesis twice in
#'   observational studies. Biometrika 99(4), 763-774. See also the CRAN package
#'   `testtwice`.
#' @export
riposte_pairs_bound <- function(y, H, gamma = 1, algorithm = NULL) {
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("package `mvtnorm` is required.", call. = FALSE)
  H <- as.matrix(H)
  if (nrow(H) != length(y)) stop("`H` needs one row per element of `y`.", call. = FALSE)
  if (ncol(H) < 2L && !isTRUE(attr(H, "riposte_subset")))
    stop("`H` needs at least two columns; use a single statistic directly.", call. = FALSE)
  if (!is.numeric(gamma) || length(gamma) != 1L || gamma < 1) stop("`gamma` must be a single number >= 1.", call. = FALSE)
  if (any(H < 0)) stop("scores must be nonnegative.", call. = FALSE)
  k <- ncol(H)
  nm <- colnames(H) %||% paste0("T", seq_len(k))
  ## a score column that is identically zero carries no information and has no variance:
  ## leave it out before standardizing; if every column is zero there is no evidence at all
  live <- colSums(H) > 0
  if (!all(live)) {
    dev <- stats::setNames(rep(NA_real_, k), nm)
    if (!any(live))
      return(list(p.value = 1, deviates = dev, cor = matrix(NA_real_, k, k, dimnames = list(nm, nm))))
    Hl <- H[, live, drop = FALSE]; attr(Hl, "riposte_subset") <- TRUE
    sub <- riposte_pairs_bound(y, Hl, gamma = gamma, algorithm = algorithm)
    dev[live] <- sub$deviates
    cr <- matrix(NA_real_, k, k, dimnames = list(nm, nm))
    cr[live, live] <- sub$cor
    return(list(p.value = sub$p.value, deviates = dev, cor = cr))
  }
  B <- as.numeric(y > 0)
  theta <- gamma / (1 + gamma)
  st <- as.vector(B %*% H)
  ex <- theta * colSums(H)
  cv <- theta * (1 - theta) * crossprod(H)
  dev <- (st - ex) / sqrt(diag(cv))
  cr <- stats::cov2cor(cv)
  # a statistic that is a perfect copy of an earlier one (this happens with few pairs)
  # adds nothing to the maximum but makes the correlation matrix singular, which Miwa
  # cannot handle; leave the copies out of the probability (the maximum itself is unchanged)
  dup <- vapply(seq_len(k), function(j) any(abs(cr[j, seq_len(j - 1L)] - 1) < 1e-10),
                logical(1))
  crk <- cr[!dup, !dup, drop = FALSE]
  kk <- ncol(crk)
  if (is.null(algorithm))
    algorithm <- if (kk <= 5 && kk > 1 &&
                     min(eigen(crk, symmetric = TRUE, only.values = TRUE)$values) > 1e-8)
      mvtnorm::Miwa(steps = 512) else mvtnorm::GenzBretz()
  p <- if (kk == 1L) 1 - stats::pnorm(max(dev)) else
    1 - mvtnorm::pmvnorm(lower = rep(-Inf, kk), upper = rep(max(dev), kk),
                         corr = crk, algorithm = algorithm)[[1]]
  names(dev) <- nm; dimnames(cr) <- list(nm, nm)
  list(p.value = p, deviates = dev, cor = cr)
}

#' Test one hypothesis with several matched-pair statistics
#'
#' Convenience wrapper: builds the score columns with [riposte_pairs_scores()] and
#' calls [riposte_pairs_bound()].
#'
#' @param y matched-pair differences (treated minus control).
#' @param statistics character vector from `"wilcoxon"`, `"brown"`, `"noether"`,
#'   `"u858"`, `"u888"`, `"u878"`, `"u868"`, `"u867"` (U-statistics with `m = 8`).
#'   Default: the two U-statistics `u878` and `u867` highlighted in Rosenbaum (2012).
#' @param gamma sensitivity parameter.
#' @param alternative `"greater"` or `"less"`.
#' @param dose optional nonnegative weights (e.g. doses) multiplying every score.
#' @param algorithm passed to [riposte_pairs_bound()].
#' @return a [riposte_pairs_bound()] list.
#' @export
riposte_pairs_test <- function(y, statistics = c("u878", "u867"), gamma = 1,
                               alternative = c("greater", "less"), dose = NULL,
                               algorithm = NULL) {
  alternative <- match.arg(alternative)
  if (alternative == "less") y <- -y
  if (!is.null(dose) && (length(dose) != length(y) || any(dose < 0)))
    stop("`dose` must be nonnegative with one value per pair.", call. = FALSE)
  uspec <- list(u858 = c(8, 5, 8), u888 = c(8, 8, 8), u878 = c(8, 7, 8),
                u868 = c(8, 6, 8), u867 = c(8, 6, 7))
  cols <- lapply(statistics, function(s) {
    if (s %in% names(uspec)) {
      u <- uspec[[s]]
      riposte_pairs_scores(y, "ustat", m = u[1], m1 = u[2], m2 = u[3])
    } else riposte_pairs_scores(y, match.arg(s, c("wilcoxon", "brown", "noether")))
  })
  H <- do.call(cbind, cols); colnames(H) <- statistics
  if (!is.null(dose)) H <- H * dose
  riposte_pairs_bound(y, H, gamma = gamma, algorithm = algorithm)
}

#' Sign-pattern test for the maximum of L statistics, with a sensitivity bound
#'
#' Under the sharp null in a paired randomized experiment, each difference's sign is
#' equally likely to be either, independently across pairs. Under Rosenbaum's
#' sensitivity model with hidden bias at most `gamma`, each sign is instead positive
#' with some chance no larger than `theta = gamma / (1 + gamma)`. Because the scores
#' are nonnegative, the chance that the maximum standardized statistic reaches its
#' observed value is largest when every sign is positive with chance `theta`, so that
#' chance is an upper bound on the p-value (and the exact p-value at `gamma = 1`).
#' This function computes it from the distribution of the signs: by listing every
#' sign pattern with its chance when `nresample` is `NULL` (at most 16 pairs with
#' nonzero difference), otherwise from `nresample` random patterns with the observed
#' signs counted as one of them. [riposte_pairs_bound()] gives the normal
#' approximation to the same bound.
#'
#' @inheritParams riposte_pairs_bound
#' @param nresample number of random sign vectors, or `NULL` for full enumeration
#'   (requires at most 16 nonzero pairs).
#' @param gamma sensitivity parameter, `>= 1`; `1` (the default) gives the
#'   randomization test of the sharp null.
#' @return a list with `p.value` (an upper bound when `gamma > 1`) and the observed
#'   `statistic`, the maximum of the statistics standardized at `theta`.
#' @references Rosenbaum, P. R. (2012) Testing one hypothesis twice in
#'   observational studies. Biometrika 99(4), 763-774.
#' @export
riposte_pairs_perm_test <- function(y, H, nresample = 9999L, gamma = 1) {
  riposte_check_gamma(gamma)
  theta <- gamma / (1 + gamma)
  H <- as.matrix(H)
  keep <- y != 0
  y <- y[keep]; H <- H[keep, , drop = FALSE]
  I <- length(y)
  ## drop all-zero score columns (no variance); with none left there is no evidence
  H <- H[, colSums(H^2) > 0, drop = FALSE]
  if (ncol(H) == 0L) return(list(p.value = 1, statistic = NA_real_))
  ## mean and sd of each statistic when every sign is positive with chance theta
  mu <- theta * colSums(H)
  sdv <- sqrt(theta * (1 - theta) * colSums(H^2))
  stat_of <- function(S) {                   # S: I x R matrix of 0/1 signs
    T <- crossprod(S, H)                     # R x L
    Z <- sweep(sweep(T, 2, mu), 2, sdv, "/")
    apply(Z, 1, max)
  }
  obs <- stat_of(matrix(as.numeric(y > 0), ncol = 1))
  if (is.null(nresample)) {
    if (I > 16) stop("full enumeration needs <= 16 nonzero pairs; set `nresample`.", call. = FALSE)
    S <- as.matrix(expand.grid(rep(list(0:1), I)))
    S <- t(S)
    sims <- stat_of(S)
    ## each pattern's chance: theta for every positive sign, 1 - theta for every negative
    npos <- colSums(S)
    chance <- theta^npos * (1 - theta)^(I - npos)
    p <- sum(chance[sims >= obs - 1e-12])
  } else {
    S <- matrix(stats::rbinom(I * nresample, 1, theta), nrow = I)
    sims <- stat_of(S)
    p <- (1 + sum(sims >= obs - 1e-12)) / (nresample + 1)
  }
  list(p.value = p, statistic = unname(obs))
}

## The exact adaptive test (Rosenbaum 2012, Ann. Appl. Stat.)

#' Probabilities behind the exact adaptive test
#'
#' For band sizes `I1`, `I2` and `kappa = gamma / (1 + gamma)`, the exact joint
#' distribution of `B1` and `T = 2*B1 + B2` under the sensitivity bound gives, for
#' every pair of cutoffs, the chance of rejecting.
#'
#' @return a list: `U[kB + 1, kT + 1]` = P(B1 >= kB or T >= kT); `pB[kB + 1]` =
#'   P(B1 >= kB); `pT[kT + 1]` = P(T >= kT).
#' @keywords internal
#' @noRd
riposte_adaptive_grid <- function(I1, I2, gamma = 1) {
  riposte_check_gamma(gamma)
  kappa <- gamma / (1 + gamma)
  p1 <- stats::dbinom(0:I1, I1, kappa)       # P(B1 = 0), ..., P(B1 = I1)
  cdf2 <- function(x) ifelse(x < 0, 0, stats::pbinom(pmin(x, I2), I2, kappa))
  kTmax <- 2 * I1 + I2 + 1L                  # one above the largest possible T
  none <- matrix(0, I1 + 2L, kTmax + 1L)     # P(B1 < kB and T < kT)
  for (kT in 0:kTmax)
    none[, kT + 1L] <- c(0, cumsum(p1 * cdf2(kT - 2 * (0:I1) - 1)))
  list(U = 1 - none, pB = 1 - c(0, cumsum(p1)), pT = 1 - none[I1 + 2L, ],
       I1 = I1, I2 = I2, kTmax = kTmax)
}

#' Choose the cutoffs (kB, kT) at level alpha from the precomputed grid
#' @keywords internal
#' @noRd
riposte_adaptive_pick <- function(g, alpha) {
  I1 <- g$I1; kTmax <- g$kTmax
  kB <- 1:(I1 + 1L); kT <- 1:kTmax
  here  <- g$U[kB + 1L, kT + 1L, drop = FALSE]     # size at (kB, kT)
  lowB  <- g$U[kB, kT + 1L, drop = FALSE]          # size at (kB - 1, kT)
  lowT  <- g$U[kB + 1L, kT, drop = FALSE]          # size at (kB, kT - 1)
  ok <- here <= alpha & lowB > alpha & lowT > alpha         # conditions (4) and (5)
  if (!any(ok)) return(NULL)
  d <- abs(outer(g$pB[kB + 1L], g$pT[kT + 1L], "-"))        # condition (6)
  d[!ok] <- Inf
  dmin <- min(d)
  ## exact ties go to the smallest kB, as in the paper
  hit <- which(d <= dmin + 1e-15, arr.ind = TRUE)
  hit <- hit[order(hit[, 1], hit[, 2]), , drop = FALSE][1, ]
  list(kB = kB[hit[1]], kT = kT[hit[2]], size = here[hit[1], hit[2]],
       pB = g$pB[kB[hit[1]] + 1L], pT = g$pT[kT[hit[2]] + 1L])
}

#' Critical values of the exact adaptive test
#'
#' With `I1` pairs in the top `lambda[1]` fraction of ranks of `|y|` and `I2` in the
#' next band, `B1`, `B2` count positive differences, and `T = 2*B1 + B2`. Under the
#' sensitivity bound they are independent binomials with success probability
#' `kappa = gamma / (1 + gamma)`. The test rejects when `B1 >= kB` or `T >= kT`,
#' where `(kB, kT)` satisfy conditions (4)-(6) of the paper: size at most `alpha`,
#' neither constant can be lowered by one without exceeding `alpha`, and
#' `|P(B1 >= kB) - P(T >= kT)|` is minimized. Ties in the minimization go to the
#' smaller `kB`, as in the paper.
#'
#' @param I1,I2 band sizes.
#' @param gamma sensitivity parameter, `>= 1`.
#' @param alpha level.
#' @return a list with `kB`, `kT`, `size` (P(reject)), `pB` (P(B1 >= kB)), `pT`.
#' @export
riposte_adaptive_critical <- function(I1, I2, gamma = 1, alpha = 0.05) {
  riposte_check_gamma(gamma)
  riposte_check_alpha(alpha)
  cr <- riposte_adaptive_pick(riposte_adaptive_grid(I1, I2, gamma), alpha)
  if (is.null(cr)) stop("no (kB, kT) satisfies the conditions; sample too small for this alpha/gamma.", call. = FALSE)
  cr
}

#' Split pairs into the adaptive test's two rank bands and count positives
#' @keywords internal
#' @noRd
riposte_adaptive_counts <- function(y, lambda = c(1 / 3, 2 / 3)) {
  riposte_check_lambda(lambda)
  y <- y[y != 0]
  I <- length(y)
  if (anyDuplicated(abs(y))) warning("ties among |y|: broken by order; the paper assumes untied ranks.", call. = FALSE)
  q <- rank(abs(y), ties.method = "first")
  in1 <- riposte_rank_at_least(q, I, 1 - lambda[1])
  in2 <- !in1 & riposte_rank_at_least(q, I, 1 - lambda[2])
  B1 <- sum(y[in1] > 0); B2 <- sum(y[in2] > 0)
  list(I1 = sum(in1), I2 = sum(in2), B1 = B1, T = 2 * B1 + B2)
}

#' P-value of the exact adaptive test
#'
#' Write `alpha*` for the smallest level at which the adaptive test rejects the
#' observed counts. Rosenbaum chooses the cutoffs `(kB, kT)` separately at each
#' level, so the rejection regions at different levels are not always nested:
#' counts rejected at one level can go unrejected at a larger one, and `alpha*`
#' itself, reported as a p-value, would fall at or below a level more often than
#' that level allows. The reported p-value treats `alpha*` as a test statistic. It
#' is the probability, under the sensitivity bound, of counts whose smallest
#' rejecting level is at most the observed `alpha*`, which is the probability of the
#' union of the rejection regions at every level up to `alpha*`. It is exact at
#' `gamma = 1` and an upper bound for `gamma > 1`, so `P(p.value <= a) <= a`. When
#' the regions are nested it is at most `alpha*`.
#'
#' The chance of each rejection region is computed as 1 minus the chance of the
#' counts outside it, and that subtraction loses all its digits once the region's
#' chance is near `1e-16`. Below about `1e-15`, rounding error rather than
#' Rosenbaum's conditions decides the cutoffs. A reported p-value that small is
#' still valid, because it is the exact chance of the region the code chose, but it
#' should be read only as "below about `1e-15`".
#'
#' @inheritParams riposte_adaptive_test
#' @param scan not used; kept so that calls written for an earlier version still run.
#' @return a single p-value in (0, 1].
#' @export
riposte_adaptive_pvalue <- function(y, gamma = 1, lambda = c(1 / 3, 2 / 3), scan = NULL) {
  riposte_check_gamma(gamma)
  ct <- riposte_adaptive_counts(y, lambda)
  g <- riposte_adaptive_grid(ct$I1, ct$I2, gamma)
  kappa <- gamma / (1 + gamma)
  ## null probability of every (B1, B2) cell, and which cells each cutoff pair rejects
  pr <- outer(stats::dbinom(0:ct$I1, ct$I1, kappa), stats::dbinom(0:ct$I2, ct$I2, kappa))
  Tm <- outer(2 * (0:ct$I1), 0:ct$I2, "+")
  B1m <- matrix(0:ct$I1, ct$I1 + 1L, ct$I2 + 1L)
  ## candidate levels: every achievable size of a rejection region (all of them, not
  ## just those containing the observed counts, because the union must include the
  ## regions of every smaller level too)
  kB <- 1:(ct$I1 + 1L); kT <- 1:g$kTmax
  cand <- sort(unique(g$U[kB + 1L, kT + 1L]))
  cand <- cand[cand <= 1]
  union <- matrix(FALSE, ct$I1 + 1L, ct$I2 + 1L)
  for (a in cand) {
    cr <- riposte_adaptive_pick(g, a)
    if (is.null(cr)) next
    union <- union | B1m >= cr$kB | Tm >= cr$kT
    if (ct$B1 >= cr$kB || ct$T >= cr$kT) return(min(1, sum(pr[union])))
  }
  1
}

#' Exact adaptive test for matched pairs (Brown + Noether)
#'
#' @param y matched-pair differences (treated minus control). Zero differences are
#'   dropped. Ties among `|y|` are broken by order (`ties.method = "first"`) with a
#'   warning, because the paper assumes untied ranks.
#' @param gamma sensitivity parameter.
#' @param alpha level.
#' @param lambda the two band fractions `0 < lambda[1] < lambda[2] <= 1`
#'   (default `c(1/3, 2/3)`).
#' @return a list with `reject`, `B1`, `T`, `kB`, `kT`, `I1`, `I2`, `size`.
#' @references Rosenbaum, P. R. (2012) An exact adaptive test with superior design
#'   sensitivity in an observational study of treatments for ovarian cancer. Annals
#'   of Applied Statistics 6(1), 83-105.
#' @export
riposte_adaptive_test <- function(y, gamma = 1, alpha = 0.05, lambda = c(1 / 3, 2 / 3)) {
  ct <- riposte_adaptive_counts(y, lambda)             # bands and positive counts
  cr <- riposte_adaptive_critical(ct$I1, ct$I2, gamma, alpha)   # the cutoffs kB and kT
  list(reject = (ct$B1 >= cr$kB) || (ct$T >= cr$kT), B1 = ct$B1, T = ct$T,
       kB = cr$kB, kT = cr$kT, I1 = ct$I1, I2 = ct$I2, size = cr$size)
}
