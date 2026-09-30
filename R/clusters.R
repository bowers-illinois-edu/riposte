## clusters.R --- cluster-randomized designs by collapsing to the cluster level.
##
## THE POINT. Under cluster randomization, treatment is assigned to whole
## clusters, so the randomization is over CLUSTER-level assignments and the
## effective sample size is the number of clusters, not units. A cluster-
## randomized test of the sharp null is therefore exactly a block-randomized test
## whose sampling units are the clusters: collapse each cluster to a summary
## outcome (its mean by default) and summary covariates (cluster means), map each
## cluster to its block, and run the same machinery riposte uses for blocks. The
## permutation then draws cluster assignments within blocks, the closed-form
## moments are exact at the cluster level, and covariance-adjustment folds fall at
## the cluster level automatically (the learner fits on control CLUSTERS).
##
## This reduction is exact for the sharp null at ANY number of clusters -- the
## permutation reference is the true cluster randomization. (The few-cluster
## variance for ATE/Neyman estimation is a separate, open problem handled by
## propertee, not by this Fisher test.)
##
## Two design facts the collapse must enforce: treatment is constant within a
## cluster (that is what "cluster-randomized" means), and clusters nest within
## blocks (a cluster cannot straddle two blocks). Both are validated, with a clear
## error naming the offending cluster.

#' Resolve a cluster specification to a vector of cluster labels
#'
#' @param clusters a vector of cluster labels (length matching the data) or a
#'   single column name in `data`.
#' @param data the data frame.
#' @param n the number of units (for a length check).
#' @return a vector of cluster labels, length `n`.
#' @keywords internal
#' @noRd
riposte_resolve_clusters <- function(clusters, data, n) {
  cl <- if (length(clusters) == 1L && is.character(clusters) &&
            clusters %in% names(data)) data[[clusters]] else clusters
  if (length(cl) != n)
    stop("cluster length does not match the data.", call. = FALSE)
  if (anyNA(cl))
    stop(sprintf("cluster has %d missing value(s); remove or impute them before testing.",
                 sum(is.na(cl))), call. = FALSE)
  cl
}

#' Collapse a unit-level design to the cluster level
#'
#' Aggregates the outcome (by `agg`, the cluster mean by default) and the
#' covariates (cluster means) to one row per cluster, carrying each cluster's
#' (constant) treatment and (unique) block. Validates that treatment is constant
#' within each cluster and that clusters nest within blocks.
#'
#' @param y,z,block unit-level outcome, 0/1 treatment, and block factor.
#' @param cluster unit-level cluster labels.
#' @param X optional unit-level covariate matrix (collapsed by cluster mean).
#' @param agg outcome aggregation function (default [mean()]).
#' @return a list with cluster-level `y`, `z`, `block`, `X` (or `NULL`), the
#'   `cluster` factor, `n_clusters`, and `n_units`.
#' @keywords internal
#' @noRd
riposte_collapse_clusters <- function(y, z, block, cluster, X = NULL, agg = mean) {
  cluster <- as.factor(cluster)
  idx <- split(seq_along(cluster), cluster, drop = TRUE)

  ## validate the two cluster-randomization facts, naming the offending cluster
  for (g in names(idx)) {
    i <- idx[[g]]
    if (length(unique(z[i])) > 1L)
      stop(sprintf("treatment varies within cluster '%s'; cluster-randomized treatment must be constant within a cluster.",
                   g), call. = FALSE)
    if (length(unique(as.character(block[i]))) > 1L)
      stop(sprintf("cluster '%s' spans more than one block; clusters must nest within blocks.",
                   g), call. = FALSE)
  }

  yg <- vapply(idx, function(i) agg(y[i]), numeric(1))
  zg <- vapply(idx, function(i) z[i][1], numeric(1))
  bg <- factor(vapply(idx, function(i) as.character(block[i][1]), character(1)))
  Xg <- if (!is.null(X))
    do.call(rbind, lapply(idx, function(i) colMeans(X[i, , drop = FALSE]))) else NULL

  list(y = unname(yg), z = unname(zg), block = bg, X = Xg,
       cluster = factor(names(idx)), n_clusters = length(idx), n_units = length(cluster))
}
