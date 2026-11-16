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
#' @param block scores whose block statistic joins the min-p, one more p-value
#'   each. The block statistic standardizes a score's treated sum within each
#'   block, squares it, and adds the squares over the blocks, so a treatment that
#'   raises outcomes in some blocks and lowers them in others does not cancel.
#'   Give the names of representations in `representations` (for example
#'   `block = "rank"`) or a named list of score functions (for example
#'   `block = list(cosine = riposte_cosine)`, the paper's cosine statistic).
#'   `NULL` (default) adds none. Its p-value follows `pvalue`; under
#'   `"asymptotic"` it is the chi-square tail on the number of blocks in the
#'   sum. It has no one-sided form, is not yet available with covariance
#'   adjustment, and is left out, with a message, when the design has one block
#'   or only blocks of two units, where it cannot add anything.
#' @return a list of class `riposte_minp_control` with elements `pvalue`,
#'   `quadratic`, and `block`, for the `minp` argument of [riposte_test()].
#' @examples
#' riposte_minp_control()
#' riposte_minp_control(pvalue = "randomization")
#' riposte_minp_control(quadratic = FALSE)
#' riposte_minp_control(block = list(cosine = riposte_cosine))
#' @export
riposte_minp_control <- function(pvalue = c("midp", "randomization", "asymptotic"),
                                 quadratic = NULL, block = NULL) {
  pvalue <- match.arg(pvalue)
  if (!is.null(quadratic) &&
      !(is.logical(quadratic) && length(quadratic) == 1L && !is.na(quadratic)))
    stop("`quadratic` must be NULL, TRUE, or FALSE.", call. = FALSE)
  if (!is.null(block)) {
    if (is.character(block)) {
      if (length(block) == 0L || anyNA(block) || any(block == ""))
        stop("`block` must name at least one representation.", call. = FALSE)
    } else if (is.list(block)) {
      if (length(block) == 0L || is.null(names(block)) || any(names(block) == "") ||
          anyDuplicated(names(block)))
        stop("`block` must be a named list of score functions, each with its own name.",
             call. = FALSE)
      if (!all(vapply(block, is.function, logical(1))))
        stop("every element of `block` must be a function of one block's outcomes.",
             call. = FALSE)
    } else {
      stop("`block` must be NULL, the names of representations, or a named list ",
           "of score functions.", call. = FALSE)
    }
  }
  structure(list(pvalue = pvalue, quadratic = quadratic, block = block),
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

## The same error for a block statistic, which also squares (Jake, 2026-11-16).
riposte_minp_block_one_sided_error <- function(alternative) {
  moved <- if (alternative == "greater") "fell" else "rose"
  stop("The min-p cannot include a block statistic when alternative = \"",
       alternative, "\". A block statistic squares each block's standardized ",
       "comparison, so it responds to departures in both directions. A min-p that ",
       "included it could reject because the treated sums ", moved, " in some blocks. ",
       "A test of \"", alternative, "\" should not reject for that reason. Leave out ",
       "`block` in riposte_minp_control() for a one-sided min-p, or use ",
       "alternative = \"two.sided\" to include the block statistic.", call. = FALSE)
}

## the scores whose block statistics a min-p result combined, NULL for none
riposte_minp_block_names <- function(component_p) {
  nm <- grep("^block_", names(component_p), value = TRUE)
  if (length(nm)) sub("^block_", "", nm) else NULL
}

## ---- the block statistic --------------------------------------------------------

## In each block, the treated sum of the score's departures from the block mean,
## squared and divided by its variance under the design, added over the blocks
## whose variance is positive (Section 4.7 of the paper; the paper's
## member_square_then_sum()). `draws` holds one assignment per column, the
## observed one first, from which each block's treated count is read. Returns
## the statistic at every assignment and the number of blocks in the sum, the
## degrees of freedom of the chi-square reference.
riposte_block_squared <- function(score, block, draws) {
  stat <- numeric(ncol(draws))
  used <- 0L
  for (ix in split(seq_along(score), block)) {
    n_b <- length(ix)
    m_b <- sum(draws[ix, 1])
    if (n_b < 2L || m_b == 0 || m_b == n_b) next
    dev <- score[ix] - mean(score[ix])
    v_b <- m_b * (n_b - m_b) / (n_b - 1) * sum(dev^2) / n_b
    if (v_b <= 0) next
    used <- used + 1L
    stat <- stat + drop(crossprod(draws[ix, , drop = FALSE], dev))^2 / v_b
  }
  list(stat = stat, df = used)
}

## The score matrix for the block statistics that riposte_minp_control(block = )
## asks for, or NULL when none can add anything. With one block that has both
## treated and control units, the block statistic is the square of the score's
## standardized statistic: a two-sided test of the score, not a statistic that
## keeps blocks from cancelling. In blocks of two units every block's term is 1
## at every assignment, so the statistic never changes. Either way riposte
## leaves it out with a message rather than raise the number of p-values the
## min-p combines (PLAN_block_statistic_2026-11-16.md, Jake's approval).
riposte_minp_block_scores <- function(y, z, block, representations, spec) {
  funs <- if (is.character(spec)) representations[spec] else spec
  by_block <- split(z, block)
  usable <- vapply(by_block, function(v) sum(v) > 0 && sum(v) < length(v), logical(1))
  if (sum(usable) <= 1L) {
    message("The block statistic is left out: with one block it is the square of ",
            "the score's standardized statistic, a two-sided test of the score rather ",
            "than a statistic that keeps blocks from cancelling. To test the score, ",
            "add it to `representations`.")
    return(NULL)
  }
  if (all(lengths(by_block)[usable] == 2L)) {
    message("The block statistic is left out: in blocks of two units every block's ",
            "term equals 1 at every re-randomization, so the statistic never changes ",
            "and cannot reject.")
    return(NULL)
  }
  sm <- riposte_score_matrix(y, block, funs)
  if (length(sm$dropped))
    message("The block statistic is left out for ", paste(sm$dropped, collapse = ", "),
            ": the score is constant within every block.")
  if (ncol(sm$scores) == 0L) return(NULL)
  sm$scores
}

## ---- the combination on a statistics matrix ----------------------------------

## Tmat is (1 + nresample) x representations, row 1 the observed assignment. mu
## and Sigma are the null mean and covariance of the treated sums: closed-form
## without adjustment, pooled over the rows with it. The K p-values are kept on
## the log scale, so that a formula p-value too small to represent still orders
## the assignments, and the smallest is compared across rows with riposte's
## relative tie tolerance.
riposte_minp_from_T <- function(Tmat, alternative = "two.sided", pvalue = "midp",
                                quadratic = FALSE, mu = NULL, Sigma = NULL,
                                block_stats = NULL) {
  if (quadratic && alternative != "two.sided")
    riposte_minp_one_sided_error(alternative)
  if (length(block_stats) && alternative != "two.sided")
    riposte_minp_block_one_sided_error(alternative)
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
  ## each block statistic: large values are extreme
  for (nm in names(block_stats)) {
    bs <- block_stats[[nm]]
    logb <- switch(
      pvalue,
      midp = log(riposte_midp(bs$stat, alternative = "greater")),
      randomization = log(riposte_randp(bs$stat, alternative = "greater")),
      asymptotic = stats::pchisq(bs$stat, df = bs$df, lower.tail = FALSE, log.p = TRUE))
    logP <- cbind(logP, logb)
    colnames(logP)[ncol(logP)] <- paste0("block_", nm)
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
#' @param block_scores optional matrix of scores, one named column per score,
#'   whose block statistics join the min-p.
#' @return a list with the observed `statistic` (the smallest p-value), its
#'   permutation `p.value`, `which` (the p-value that was smallest),
#'   `component_p` (each observed p-value, block statistics named `block_`
#'   followed by the score's name), `pvalue`, `quadratic`, and `df` (the
#'   quadratic's rank, `NA` without it).
#' @keywords internal
#' @noRd
riposte_minp <- function(scores, z, block, moments = NULL, draws = NULL,
                         nresample = 1999L, alternative = "two.sided",
                         pvalue = "midp", quadratic = FALSE, block_scores = NULL) {
  scores <- as.matrix(scores)
  if (ncol(scores) == 0L)
    return(list(statistic = NA_real_, p.value = NA_real_, which = NA_character_,
                component_p = numeric(0)))
  moments <- moments %||% riposte_sw_moments(scores, z, block)
  perm <- riposte_perm_stats(scores, z, block, nresample, draws)
  block_stats <- NULL
  if (!is.null(block_scores)) {
    block_scores <- as.matrix(block_scores)
    block_stats <- lapply(stats::setNames(seq_len(ncol(block_scores)), colnames(block_scores)),
                          function(j) riposte_block_squared(block_scores[, j], block, perm$draws))
  }
  riposte_minp_from_T(perm$stats, alternative, pvalue, quadratic, moments$mu,
                      moments$Sigma, block_stats)
}
