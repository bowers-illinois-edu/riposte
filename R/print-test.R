## Report which reference was used, so an approximate p-value cannot be
## mistaken for a count from randomized assignments.

#' @export
print.riposte_test <- function(x, ...) {
  cat("riposte test of the sharp null of no effect\n")
  if (isTRUE(x$engine == "asymptotic")) {
    if (isTRUE(x$clustered))
      cat(sprintf("  %d units in %d clusters in %d blocks\n",
                  x$n_units, x$n_clusters, x$nblocks))
    else
      cat(sprintf("  %d units in %d blocks\n", x$n, x$nblocks))
    reference <- if (x$combination == "quadratic")
      sprintf("chi-square, df = %d", x$df) else "standard Cauchy"
    cat(sprintf("  asymptotic engine (%s; no re-randomization)\n", reference))
  } else if (isTRUE(x$engine == "saddlepoint"))
    cat(sprintf("  %d units in %d blocks; saddlepoint engine (no re-randomization)\n",
                x$n, x$nblocks))
  else if (isTRUE(x$clustered))
    cat(sprintf("  %d units in %d clusters in %d blocks; %d cluster re-randomizations\n",
                x$n_units, x$n_clusters, x$nblocks, x$nresample))
  else
    cat(sprintf("  %d units in %d blocks; %d re-randomizations\n",
                x$n, x$nblocks, x$nresample))
  combo <- x$combination
  if (x$requested == "screen") {
    how <- if (!is.null(x$screen) && x$screen$method == "shrink")
      sprintf("shrink, lambda = %.3f", x$screen$lambda)
    else sprintf("threshold, chose %s", combo)
    cat(sprintf("  screen (%s); condition number %.1f\n", how, x$condition))
  }
  cat(sprintf("  combination: %s\n", combo))
  if (isTRUE(x$adjusted))
    cat(sprintf("  covariance-adjusted (%s)\n",
                if (isTRUE(x$refit)) "refit per permutation, exact"
                else "fit once, fixed residuals -- not exact if the learner overfits"))
  if (isTRUE(x$engine == "asymptotic"))
    cat(sprintf("  statistic = %.4f,  approximate p-value = %.6g\n",
                x$statistic, x$p.value))
  else
    cat(sprintf("  statistic = %.4f,  p-value = %.4f\n", x$statistic, x$p.value))
  if (isTRUE(x$engine == "saddlepoint"))
    cat(sprintf("  (saddlepoint %s p-value; approximate, no re-randomization)\n",
                combo))
  if (length(x$dropped))
    cat("  dropped (constant under the design): ",
        paste(x$dropped, collapse = ", "), "\n")
  invisible(x)
}
