## components.R --- the per-representation diagnostics.   [STUB]
##
## THE POINT. riposte_components() exposes the building blocks the combinations
## sit on: each representation's two-sided permutation p-value, the closed-form
## permutation covariance Sigma and its condition number, and the score matrix.
## It is the diagnostic table behind the cost-of-hedging comparison and the
## screen.

#' Per-representation p-values, covariance, and condition number
#'
#' @param formula `Y ~ treatment` for complete randomization, or
#'   `Y ~ treatment | block` for randomization within blocks (or use `blocks`).
#' @param data a data.frame/data.table.
#' @param blocks optional block specification (see [riposte_test()]).
#' @param clusters optional cluster specification (see [riposte_test()]); the
#'   design is collapsed to the cluster level when given.
#' @param representations a named list of representation functions;
#'   [riposte_reps_default()] by default.
#' @param nresample number of re-randomizations for the marginal p-values.
#' @param cluster_agg how to aggregate the outcome within a cluster; the cluster
#'   mean by default.
#' @param seed optional integer seed (L'Ecuyer-CMRG) for reproducibility.
#' @param ... reserved.
#' @return an object of class `riposte_components` with the representations' two
#'   sided mid-p permutation p-values, the closed-form covariance `Sigma`, its
#'   `condition` number, the kept/dropped representations, and the score matrix.
#' @export
riposte_components <- function(formula, data, blocks = NULL, clusters = NULL,
                              representations = riposte_reps_default(),
                              nresample = 1999L, cluster_agg = mean,
                              seed = NULL, ...) {
  nresample <- riposte_check_nresample(nresample)
  if (!is.null(seed)) { RNGkind("L'Ecuyer-CMRG"); set.seed(seed) }
  des <- riposte_parse_design(formula, data, blocks)
  n_units <- length(des$y); n_clusters <- NA_integer_
  if (!is.null(clusters)) {
    cl <- riposte_resolve_clusters(clusters, data, length(des$y))
    col <- riposte_collapse_clusters(des$y, des$z, des$block, cl, NULL, cluster_agg)
    des <- list(y = col$y, z = col$z, block = col$block)
    n_clusters <- col$n_clusters
  }
  sm <- riposte_score_matrix(des$y, des$block, representations)
  if (ncol(sm$scores) == 0L)
    stop("every representation is constant within blocks; nothing to test.",
         call. = FALSE)
  moments <- riposte_sw_moments(sm$scores, des$z, des$block)
  riposte_assert_testable(moments$Sigma)
  Tmat <- riposte_perm_stats(sm$scores, des$z, des$block, nresample)$stats
  comp_p <- riposte_midp_matrix_first(Tmat)
  scr <- riposte_screen(moments$Sigma)

  structure(list(
    component_p = comp_p,
    Sigma = moments$Sigma,
    condition = scr$condition,
    kept = sm$kept,
    dropped = sm$dropped,
    scores = sm$scores,
    n = length(des$y),
    n_units = n_units,
    n_clusters = n_clusters,
    clustered = !is.null(clusters),
    nblocks = nlevels(des$block),
    nresample = nresample
  ), class = "riposte_components")
}

## observed-row mid-p p-value of each representation, from the shared draws
riposte_midp_matrix_first <- function(Tmat) {
  P <- apply(Tmat, 2, riposte_midp)
  stats::setNames(P[1, ], colnames(Tmat))
}

#' @export
print.riposte_components <- function(x, ...) {
  cat("riposte components\n")
  if (isTRUE(x$clustered))
    cat(sprintf("  %d units in %d clusters in %d blocks; %d cluster re-randomizations\n",
                x$n_units, x$n_clusters, x$nblocks, x$nresample))
  else
    cat(sprintf("  %d units in %d blocks; %d re-randomizations\n",
                x$n, x$nblocks, x$nresample))
  cat(sprintf("  condition number of the representation correlations: %.1f\n",
              x$condition))
  if (length(x$dropped))
    cat("  dropped (constant within blocks): ", paste(x$dropped, collapse = ", "), "\n")
  cat("  representation two-sided mid-p p-values:\n")
  pv <- x$component_p
  for (nm in names(pv)) cat(sprintf("    %-16s %.4f\n", nm, pv[[nm]]))
  invisible(x)
}
