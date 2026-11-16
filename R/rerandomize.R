## rerandomize.R --- within-block re-randomization (the randomization reference).
##
## THE POINT. The Cauchy combination needs Monte-Carlo draws from the
## randomization distribution (you cannot get an exact two-sided permutation
## p-value from moments alone without a probit). We draw the re-randomizations
## the way the study assigns treatment: sample() the treated units within each
## block, preserving each block's treated count. The spec is emphatic about this
## --- plain sample(), not a random-key/rank shuffle. The two are equivalent, but
## the plain form reads like the design it imitates.
##
## The first column of the returned matrix is ALWAYS the observed assignment, so
## the observed statistic is exchangeable with the draws and the permutation
## p-value is the share of columns at least as extreme.

#' Within-block re-randomizations preserving each block's treated count
#'
#' Builds a 0/1 assignment matrix whose first column is the observed treatment
#' and whose remaining `nresample` columns are independent re-randomizations,
#' each drawn by `sample()` within block so that every block keeps its observed
#' number of treated units. This is the randomization reference for the Monte
#' Carlo combinations.
#'
#' For a cluster-randomized design, pass cluster-level vectors (`z` and `block`
#' indexed by cluster, one entry per cluster); the permutation is then over
#' cluster assignments, which is the correct reference for clusters.
#'
#' @param z integer/numeric 0/1 vector of observed treatment assignment.
#' @param block factor or vector of block labels, the same length as `z`.
#' @param nresample number of re-randomizations (columns 2..nresample+1).
#' @return a numeric matrix with `length(z)` rows and `nresample + 1` columns;
#'   column 1 is `z`. Row order matches the input.
#' @export
riposte_block_draws <- function(z, block, nresample = 1999L) {
  z <- as.numeric(z)
  draw <- riposte_draw_generator(z, block)
  cbind(z, draw(nresample), deparse.level = 0)
}

## A function of b that returns b new re-randomizations as an n x b matrix. Each
## draw places every block's treated count at random positions within the block,
## exactly as the study assigns treatment. Draws come one at a time from the
## random number stream, so b draws in one call and the same b draws in several
## calls are the same assignments in the same order.
riposte_draw_generator <- function(z, block) {
  z <- as.numeric(z)
  n <- length(z)
  if (length(block) != n)
    stop("`z` and `block` must have the same length.", call. = FALSE)
  ## index sets and treated counts per block, from the observed assignment
  idx_list <- split(seq_len(n), block)
  nt <- vapply(idx_list, function(ix) sum(z[ix]), numeric(1))
  one_draw <- function() {
    g <- numeric(n)
    for (bi in seq_along(idx_list)) {
      ix <- idx_list[[bi]]
      k <- nt[bi]
      if (k > 0L) g[ix[sample.int(length(ix), k)]] <- 1
    }
    g
  }
  function(b) {
    if (b <= 0L) return(matrix(numeric(0), nrow = n, ncol = 0))
    matrix(replicate(b, one_draw()), nrow = n, ncol = b)
  }
}

## Statistics at the observed assignment and at nresample re-randomizations,
## computed a batch at a time (Jake's task of 2026-11-06). fun(G) takes an n x b
## matrix of assignments, the observed one first in the first batch, and returns
## a matrix with one row per assignment. A batch holds at most `entries` matrix
## entries (option riposte.batch_entries, default 1e7, about 80 MB of doubles),
## so memory grows with the batch rather than with nresample: 100,000 units with
## 40,000 re-randomizations would otherwise need about 32 GB. With `draws`
## given, fun is applied to that matrix as a whole, which the tests that list
## every assignment use.
riposte_over_draws <- function(z, block, nresample, fun, draws = NULL,
                               entries = getOption("riposte.batch_entries", 1e7)) {
  if (!is.null(draws)) return(as.matrix(fun(draws)))
  z <- as.numeric(z)
  per <- max(1L, floor(entries / length(z)))
  draw <- riposte_draw_generator(z, block)
  total <- nresample + 1L
  out <- vector("list", ceiling(total / per))
  done <- 0L
  i <- 0L
  while (done < total) {
    b <- min(per, total - done)
    G <- if (done == 0L) cbind(z, draw(b - 1L), deparse.level = 0) else draw(b)
    i <- i + 1L
    out[[i]] <- as.matrix(fun(G))
    done <- done + b
  }
  do.call(rbind, out)
}
