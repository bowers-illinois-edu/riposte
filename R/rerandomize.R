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
  n <- length(z)
  if (length(block) != n)
    stop("`z` and `block` must have the same length.", call. = FALSE)
  ## index sets and treated counts per block, from the observed assignment
  idx_list <- split(seq_len(n), block)
  nt <- vapply(idx_list, function(ix) sum(z[ix]), numeric(1))

  ## one re-randomization: within each block, place the block's treated count at
  ## random positions, exactly as the study assigns treatment
  one_draw <- function() {
    g <- numeric(n)
    for (bi in seq_along(idx_list)) {
      ix <- idx_list[[bi]]
      k <- nt[bi]
      if (k > 0L) g[ix[sample.int(length(ix), k)]] <- 1
    }
    g
  }

  draws <- if (nresample > 0L) replicate(nresample, one_draw()) else
    matrix(numeric(0), nrow = n, ncol = 0)
  cbind(z, draws, deparse.level = 0)
}
