## rosenbaum-pairs.R
## Rosenbaum (2012) matched-pair tests, implemented in riposte.
##
## Two procedures for matched pairs, both built on signed statistics
## T = sum_i B_i h_i, where B_i = 1{y_i > 0} and h_i is a score that depends on the
## rank of |y_i|:
##   (a) "testing one hypothesis twice" (Rosenbaum 2012, Biometrika 99:763-774):
##       standardize L statistics, refer their MAXIMUM to a multivariate normal
##       orthant probability. Reference implementation: the CRAN package `testtwice`.
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
                p <- rk / I; (p >= a) + (p >= b)
              },
              noether = {
                a <- q1 %||% (2 / 3); (rk / I >= a) * 1
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
#' For `gamma = 1` an exact sign-flip alternative is [riposte_pairs_perm_test()].
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
  if (ncol(H) < 2L) stop("`H` needs at least two columns; use a single statistic directly.", call. = FALSE)
  if (!is.numeric(gamma) || length(gamma) != 1L || gamma < 1) stop("`gamma` must be a single number >= 1.", call. = FALSE)
  if (any(H < 0)) stop("scores must be nonnegative.", call. = FALSE)
  k <- ncol(H)
  nm <- colnames(H) %||% paste0("T", seq_len(k))
  B <- as.numeric(y > 0)
  theta <- gamma / (1 + gamma)
  st <- as.vector(B %*% H)
  ex <- theta * colSums(H)
  cv <- theta * (1 - theta) * crossprod(H)
  dev <- (st - ex) / sqrt(diag(cv))
  cr <- stats::cov2cor(cv)
  if (is.null(algorithm))
    algorithm <- if (k <= 5) mvtnorm::Miwa(steps = 512) else mvtnorm::GenzBretz()
  p <- 1 - mvtnorm::pmvnorm(lower = rep(-Inf, k), upper = rep(max(dev), k),
                            corr = cr, algorithm = algorithm)[[1]]
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

#' Exact sign-flip randomization test for the maximum of L statistics (gamma = 1)
#'
#' Under the sharp null in a paired randomized experiment, each difference's sign is
#' equally likely to be either, independently across pairs. This refers the maximum
#' standardized statistic to that distribution exactly (full enumeration when
#' `I <= 16` pairs with nonzero difference and `nresample` is `NULL`, otherwise
#' Monte Carlo with the observed assignment counted).
#'
#' @inheritParams riposte_pairs_bound
#' @param nresample number of random sign vectors, or `NULL` for full enumeration
#'   (requires at most 16 nonzero pairs).
#' @return a list with `p.value` and the observed `statistic` (max standardized deviate).
#' @export
riposte_pairs_perm_test <- function(y, H, nresample = 9999L) {
  H <- as.matrix(H)
  keep <- y != 0
  y <- y[keep]; H <- H[keep, , drop = FALSE]
  I <- length(y)
  sdv <- sqrt(colSums(H^2) / 4)            # theta = 1/2: Var = (1/4) sum h^2
  stat_of <- function(S) {                   # S: I x R matrix of 0/1 signs
    T <- crossprod(S, H)                     # R x L
    Z <- sweep(sweep(T, 2, colSums(H) / 2), 2, sdv, "/")
    apply(Z, 1, max)
  }
  obs <- stat_of(matrix(as.numeric(y > 0), ncol = 1))
  if (is.null(nresample)) {
    if (I > 16) stop("full enumeration needs <= 16 nonzero pairs; set `nresample`.", call. = FALSE)
    S <- as.matrix(expand.grid(rep(list(0:1), I)))
    S <- t(S)
    sims <- stat_of(S)
    p <- mean(sims >= obs - 1e-12)
  } else {
    S <- matrix(stats::rbinom(I * nresample, 1, 0.5), nrow = I)
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
  cr <- riposte_adaptive_pick(riposte_adaptive_grid(I1, I2, gamma), alpha)
  if (is.null(cr)) stop("no (kB, kT) satisfies the conditions; sample too small for this alpha/gamma.", call. = FALSE)
  cr
}

#' Split pairs into the adaptive test's two rank bands and count positives
#' @keywords internal
#' @noRd
riposte_adaptive_counts <- function(y, lambda = c(1 / 3, 2 / 3)) {
  y <- y[y != 0]
  I <- length(y)
  if (anyDuplicated(abs(y))) warning("ties among |y|: broken by order; the paper assumes untied ranks.", call. = FALSE)
  q <- rank(abs(y), ties.method = "first")
  in1 <- q >= (1 - lambda[1]) * I
  in2 <- q < (1 - lambda[1]) * I & q >= (1 - lambda[2]) * I
  B1 <- sum(y[in1] > 0); B2 <- sum(y[in2] > 0)
  list(I1 = sum(in1), I2 = sum(in2), B1 = B1, T = 2 * B1 + B2)
}

#' P-value of the exact adaptive test: the smallest level that rejects
#'
#' As for procedure (6) in the Biometrika paper, the reported p-value is the
#' smallest level `alpha` at which the adaptive test rejects. The achievable sizes
#' form a finite set, so we search it. Rejection is nearly always monotone in
#' `alpha`; we bisect and then confirm with a local scan over the next
#' `scan` smaller candidates (`scan = Inf` scans them all, which is exact).
#'
#' @inheritParams riposte_adaptive_test
#' @param scan how many smaller candidate levels to double-check after bisecting.
#' @return a single p-value in (0, 1].
#' @export
riposte_adaptive_pvalue <- function(y, gamma = 1, lambda = c(1 / 3, 2 / 3), scan = 50) {
  ct <- riposte_adaptive_counts(y, lambda)
  g <- riposte_adaptive_grid(ct$I1, ct$I2, gamma)
  rejects <- function(alpha) {
    cr <- riposte_adaptive_pick(g, alpha)
    !is.null(cr) && (ct$B1 >= cr$kB || ct$T >= cr$kT)
  }
  ## candidate levels: sizes of the rejection regions that contain the observed
  kB <- 1:(ct$I1 + 1L); kT <- 1:g$kTmax
  inreg <- outer(kB <= ct$B1, kT <= ct$T, "|")
  cand <- sort(unique(g$U[kB + 1L, kT + 1L][inreg]))
  cand <- cand[cand <= 1]
  if (!length(cand) || !rejects(max(cand))) return(1)
  lo <- 1L; hi <- length(cand)               # invariant: rejects(cand[hi])
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    if (rejects(cand[mid])) hi <- mid else lo <- mid + 1L
  }
  ## belt and braces: look for a smaller rejecting level the bisection skipped
  if (hi > 1L) {
    below <- seq_len(hi - 1L)
    below <- utils::tail(below, if (is.finite(scan)) scan else length(below))
    for (j in below) if (rejects(cand[j])) { hi <- j; break }
  }
  cand[hi]
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
