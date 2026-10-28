## pair-reps.R
## Rosenbaum's matched-pair statistics as riposte representations.
##
## A matched pair is a block of two units with one treated. Rosenbaum's
## statistic T = sum_i B_i h_i (B_i = 1 if the treated unit of pair i did better,
## h_i a score from the rank of |difference| among ALL pairs) is a Strasser-Weber
## linear statistic: give the better-off unit of each pair the score h_i and the
## other unit 0; summing scores over the treated units then gives exactly
## sum_i B_i h_i. So with these as representations riposte's closed-form moments,
## quadratic / Cauchy / max combinations, screen, and exact re-randomization all
## apply unchanged. The only new ingredient is that h_i looks across blocks, so the
## representation takes (y, block) -- class "riposte_pooled_rep" -- instead of one
## block's outcomes.

#' Matched-pair representation (Rosenbaum 2012)
#'
#' A representation for [riposte_test()] and [riposte_components()] that scores each
#' pair by the rank of its absolute difference among all pairs, and gives that score
#' to the unit with the larger outcome. Every block must be a pair (two units, one
#' treated). See [riposte_pairs_scores()] for the scores.
#'
#' @inheritParams riposte_pairs_scores
#' @return a function of `(y, block)` of class `"riposte_pooled_rep"`, carrying
#'   the pair-level score function as attribute `pair_scores`.
#' @export

riposte_pair_rep <- function(type = c("wilcoxon", "brown", "noether", "ustat"),
                             m = 8, m1 = 7, m2 = 8, exact = TRUE) {
  type <- match.arg(type)
  h_of <- function(d) riposte_pairs_scores(d, type, m = m, m1 = m1, m2 = m2, exact = exact)
  f <- function(y, block) {
    ix <- split(seq_along(y), block)
    if (any(lengths(ix) != 2L))
      stop("pair representations need matched pairs: every block must have exactly two units.",
           call. = FALSE)
    first  <- vapply(ix, function(i) i[1L], integer(1))
    second <- vapply(ix, function(i) i[2L], integer(1))
    d <- y[first] - y[second]                  # sign depends on the order; the scores do not
    h <- h_of(d)
    out <- numeric(length(y))
    out[first]  <- h * (d > 0)                 # the better-off unit of each pair gets h
    out[second] <- h * (d < 0)
    out
  }
  structure(f, class = c("riposte_pooled_rep", "function"), pair_scores = h_of)
}

#' Rosenbaum's matched-pair representations
#'
#' Wilcoxon, Brown, Noether and the two U-statistics `(8, 7, 8)` and `(8, 6, 7)`
#' highlighted in Rosenbaum (2012), as a named list for the `representations`
#' argument of [riposte_test()].
#'
#' @param which names to keep, from `"wilcoxon"`, `"brown"`, `"noether"`,
#'   `"u878"`, `"u867"`, `"u858"`, `"u888"`, `"u868"`.
#' @return a named list of [riposte_pair_rep()] functions.
#' @export
riposte_pair_reps <- function(which = c("wilcoxon", "brown", "noether", "u878", "u867")) {
  spec <- list(
    wilcoxon = riposte_pair_rep("wilcoxon"), brown = riposte_pair_rep("brown"),
    noether = riposte_pair_rep("noether"),
    u858 = riposte_pair_rep("ustat", 8, 5, 8), u888 = riposte_pair_rep("ustat", 8, 8, 8),
    u878 = riposte_pair_rep("ustat", 8, 7, 8), u868 = riposte_pair_rep("ustat", 8, 6, 8),
    u867 = riposte_pair_rep("ustat", 8, 6, 7))
  bad <- setdiff(which, names(spec))
  if (length(bad)) stop("unknown pair representation(s): ", paste(bad, collapse = ", "), call. = FALSE)
  spec[which]
}

#' Treated-minus-control difference in each matched pair
#' @keywords internal
#' @noRd
riposte_pair_differences <- function(des) {
  ix <- split(seq_along(des$y), des$block)
  ok <- all(lengths(ix) == 2L) && all(vapply(ix, function(i) sum(des$z[i]) == 1, logical(1)))
  if (!ok)
    stop("this needs matched pairs: every block must have exactly two units, one treated and one control.",
         call. = FALSE)
  vapply(ix, function(i) des$y[i][des$z[i] == 1] - des$y[i][des$z[i] == 0], numeric(1))
}

#' The Rosenbaum sensitivity-analysis branch of riposte_test()
#'
#' Called when `gamma` is supplied. One-sided (`alternative`), matched pairs only.
#' `statistic = "max"` is Rosenbaum's "testing one hypothesis twice" (procedure 6)
#' over the pair representations; `"adaptive"` is the exact Brown + Noether
#' adaptive test. With the permutation engine the max is referred to the
#' distribution of the signs when each is positive with chance gamma / (1 + gamma),
#' by listing or drawing sign patterns (exact for the sharp null at gamma = 1, an
#' upper bound on the p-value above it); with `engine = "asymptotic"`, to the
#' multivariate normal orthant probability.
#' @keywords internal
#' @noRd
riposte_test_pairs <- function(des, representations, statistic, gamma, alpha,
                               alternative, engine, nresample) {
  d <- riposte_pair_differences(des)
  if (alternative == "less") d <- -d
  base <- list(condition = NA_real_, screen = NULL, df = NA_integer_,
               dropped = character(0), gamma = gamma, alternative = alternative,
               pair_test = TRUE)
  if (statistic == "adaptive") {
    ad <- riposte_adaptive_test(d, gamma = gamma, alpha = alpha)
    return(c(list(statistic = ad$T, p.value = riposte_adaptive_pvalue(d, gamma),
                  combination = "adaptive", kept = "adaptive",
                  adaptive = c(ad, list(alpha = alpha))), base))
  }
  if (statistic != "max")
    stop("with `gamma`, statistic must be \"max\" (testing one hypothesis twice) or ",
         "\"adaptive\".", call. = FALSE)
  hs <- lapply(representations, function(f) {
    g <- attr(f, "pair_scores")
    if (is.null(g)) stop("with `gamma`, `representations` must come from riposte_pair_reps().", call. = FALSE)
    g(d)
  })
  H <- do.call(cbind, hs)
  if (ncol(H) < 2L) stop("testing twice needs at least two pair representations.", call. = FALSE)
  if (engine == "permute") {
    ## enumerate every sign pattern when that is no more work than the requested draws
    k <- sum(d != 0)
    enumerate <- k <= 16 && 2^k <= nresample
    pt <- riposte_pairs_perm_test(d, H, nresample = if (enumerate) NULL else nresample,
                                  gamma = gamma)
    return(c(list(statistic = pt$statistic, p.value = pt$p.value, combination = "max",
                  kept = names(representations), enumerated = enumerate), base))
  }
  ## engine = "asymptotic": the normal approximation of the paper and of testtwice
  bnd <- riposte_pairs_bound(d, H, gamma = gamma)
  c(list(statistic = if (all(is.na(bnd$deviates))) NA_real_ else max(bnd$deviates, na.rm = TRUE), p.value = bnd$p.value, combination = "max",
         kept = names(representations), deviates = bnd$deviates), base)
}

#' Entry for riposte_test() when `gamma` is given or statistic = "adaptive"
#' @keywords internal
#' @noRd
riposte_test_pair_design <- function(formula, data, blocks, clusters, statistic,
                                     representations, nresample, adjust, seed,
                                     engine, gamma, alpha, alternative, call) {
  if (!is.null(clusters) || !is.null(adjust))
    stop("the matched-pair tests (`gamma` / statistic = \"adaptive\") do not support ",
         "`clusters` or `adjust`.", call. = FALSE)
  if (engine == "saddlepoint")
    stop("the matched-pair tests do not use engine = \"saddlepoint\".", call. = FALSE)
  if (is.null(gamma)) gamma <- 1
  riposte_check_gamma(gamma)
  riposte_check_alpha(alpha)
  if (!is.null(seed)) { RNGkind("L'Ecuyer-CMRG"); set.seed(seed) }
  des <- riposte_parse_design(formula, data, blocks)
  out <- riposte_test_pairs(des, representations, statistic, gamma, alpha,
                            alternative, engine, nresample)
  structure(c(out, list(
    requested = statistic, engine = engine, n = length(des$y),
    n_units = length(des$y), n_clusters = NA_integer_,
    nblocks = nlevels(des$block), nresample = nresample,
    adjusted = FALSE, clustered = FALSE, call = call)), class = "riposte_test")
}
