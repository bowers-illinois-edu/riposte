## adjust.R --- controls-only covariance adjustment (Regression Without Regrets).
##
## THE POINT. Fit a learner on the CONTROL units only, residualize the outcome,
## and refer the combined statistic on the residuals to the randomization
## distribution by REFITTING the learner inside every re-randomization. Refitting
## per permutation is what keeps the test exact for ANY learner: under the sharp
## null the outcomes are fixed, so the statistic is a deterministic function of
## (assignment, outcomes, covariates), and its randomization distribution is
## exact regardless of how the learner behaves. The fit-once / fixed-residual
## shortcut breaks exactness when the learner overfits, because the residuals then
## carry information about the observed assignment (see the size tests).
##
## WHY THE COMBINED ADJUSTED TEST IS EXACT (resolved; see
## dev/theory-adjusted-exactness.md for the proposition, proof, and simulation
## evidence). The closed-form Strasser-Weber moments do NOT apply here --- the
## residual scores change from draw to draw, so the metric mu/Sigma must be
## estimated. But exact LEVEL does not need the closed-form metric. Under the
## sharp null Y is fixed, so refitting per draw makes T(z) a deterministic
## function of z; the observed and the draws are then exchangeable; and if the
## metric (mu, Sigma) is a SYMMETRIC function of the full set of statistics
## (observed + draws POOLED), the combined statistic is permutation-equivariant
## and its p-value is exactly valid for any nresample and any measurable learner.
## Two consequences:
##   - Use POOLED moments, not draws-only. Draws-only leaves the observed
##     out-of-sample to its own metric and breaks the symmetry -> anti-conservative
##     (size ~0.12 vs 0.05 at small nresample). riposte_adjusted_combination pools.
##   - The SCREEN stays exact under adjustment. Its condition number is no longer
##     "ancillary" in the closed-form sense (it depends on the realized draws), but
##     it is a symmetric function of the exchangeable pooled set, which is the
##     property the validity argument actually needs.
## "Exact" means exact LEVEL (finite-sample validity), not that the p-value equals
## the full-enumeration value: it is a Monte-Carlo permutation p-value, valid for
## any nresample, whose metric noise affects power but not level.
##
## CLUSTERS and the few-cluster SE are deferred (propertee handles the
## design-based SE; the cluster variance is open work).

#' Covariance-adjustment specification
#'
#' Bundles the covariates, the learner, and whether to refit per permutation.
#' Pass the result as `adjust =` to [riposte_test()], or pass a bare one-sided
#' covariate formula and let [riposte_test()] wrap it with the default learner.
#'
#' @param covariates a one-sided formula of covariates, e.g. `~ x1 + x2`.
#' @param learner a learner `(X, Y, control) -> Yhat`; ridge by default. See
#'   [riposte_ridge_learner()].
#' @param refit refit the learner inside every re-randomization (`TRUE`, exact
#'   for any learner) or fit once on the observed controls and permute the fixed
#'   residuals (`FALSE`, the shortcut that breaks exactness when the learner
#'   overfits --- kept for comparison).
#' @return a `riposte_adjust` specification.
#' @export
riposte_adjust <- function(covariates, learner = riposte_ridge_learner(), refit = TRUE) {
  if (!inherits(covariates, "formula"))
    stop("`covariates` must be a one-sided formula, e.g. ~ x1 + x2.", call. = FALSE)
  structure(list(covariates = covariates, learner = learner, refit = refit),
            class = "riposte_adjust")
}

#' Build the covariate matrix from a one-sided formula
#' @keywords internal
#' @noRd
riposte_covariate_matrix <- function(covariates, data) {
  mm <- stats::model.matrix(covariates, as.data.frame(data))
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

#' Residual representation statistics across assignments, with refit per draw
#'
#' For each assignment (column of `draws`, column 1 observed): residualize `y`
#' by the learner fit on that assignment's control units, build the
#' within-block-centred representation scores from the residuals, and form each
#' representation's linear statistic at that assignment. Returns the
#' (1 + nresample) x representation matrix the combinations consume. With
#' `refit = FALSE` the learner is fit once on the observed controls and the fixed
#' residual scores are reused for every assignment (the shortcut).
#'
#' @return a list with `stats` (the matrix, row 1 observed) and `kept` (the
#'   representation names retained, from the observed residuals).
#' @keywords internal
#' @noRd
riposte_residual_perm_stats <- function(y, X, block, draws, learner,
                                        representations = riposte_reps_default(),
                                        refit = TRUE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  nd <- ncol(draws)

  if (!refit) {
    ## fit once on the observed controls (column 1), permute the fixed residuals
    e <- y - learner(X, y, draws[, 1] == 0)
    sm <- riposte_score_matrix(e, block, representations)
    return(list(stats = crossprod(draws, sm$scores), kept = sm$kept))
  }

  ## Compute every representation's linear statistic at every assignment (column 1
  ## observed), refitting the learner per draw. We do NOT drop columns based on the
  ## observed assignment alone: which representations to keep must be a function of
  ## the whole orbit, not of the observed draw, or the column set -- and hence the
  ## combined statistic -- would depend on which assignment is "observed" and the
  ## exactness argument (exchangeability of the B+1 statistics) would break.
  nm <- names(representations)
  stats <- matrix(0, nd, length(representations), dimnames = list(NULL, nm))
  for (b in seq_len(nd)) {
    zb <- draws[, b]
    e <- y - learner(X, y, zb == 0)
    Sb <- riposte_centred_scores(e, block, representations)
    stats[b, ] <- colSums(Sb[zb == 1, , drop = FALSE])
  }

  ## symmetric selection: drop a representation only if its linear statistic is
  ## constant across the pooled set of assignments (zero pooled variance) -- it
  ## then carries no information and would be a degenerate column. This selection
  ## is a function of all B+1 assignments, so it preserves the exchangeability the
  ## exactness argument needs. It also leaves the pooled covariance with no
  ## zero-variance column for the screen's cov2cor.
  keep <- apply(stats, 2, function(col) stats::var(col) > 1e-12)
  if (!any(keep))
    stop("every representation is constant across the residual re-randomizations; nothing to test.",
         call. = FALSE)
  list(stats = stats[, keep, drop = FALSE], kept = nm[keep])
}

#' Run a combination on residual statistics, with Monte-Carlo moments
#'
#' @keywords internal
#' @noRd
riposte_adjusted_combination <- function(Tmat, statistic, screen_control,
                                         alternative = "two.sided") {
  ## POOLED moments (observed + draws). mu and Sigma are symmetric functions of
  ## the full set of statistics, so the combined statistic is permutation-
  ## equivariant and the test stays EXACT under refit-per-permutation (see the
  ## header and dev/theory-adjusted-exactness.md). Using the draws only would make
  ## the observed statistic out-of-sample to its own metric and the test
  ## anti-conservative (badly so at small nresample: size ~0.12 vs 0.05).
  mu <- colMeans(Tmat)
  Sigma <- stats::cov(Tmat)
  riposte_assert_testable(Sigma)

  scr <- riposte_screen(Sigma, screen_control)
  chosen <- statistic
  if (statistic == "screen")
    chosen <- if (scr$method == "threshold") scr$choice else "quadratic"

  res <- switch(
    chosen,
    quadratic = {
      Sig <- if (statistic == "screen" && scr$method == "shrink") scr$Sigma else Sigma
      riposte_quadratic_from_T(Tmat, mu, Sig)
    },
    cauchy = riposte_cauchy_from_T(Tmat, alternative),
    max    = riposte_max_from_T(Tmat, mu, sqrt(diag(Sigma)), alternative)
  )
  list(result = res, chosen = chosen, screen = if (statistic == "screen") scr else NULL,
       condition = scr$condition)
}
