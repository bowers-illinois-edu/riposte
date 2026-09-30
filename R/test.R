## test.R --- the main user entry point.
##
## THE POINT. riposte_test() runs one combined test of the sharp null of no
## effect. The default p-value uses re-randomization. statistic = "screen" (the
## default) uses the closed-form Sigma to choose/blend the quadratic and the
## Cauchy; "quadratic", "cauchy", and "max" force a single combination. adjust /
## learner add controls-only covariance adjustment. Works for complete,
## block- and cluster-randomized designs (cluster: permutation over cluster
## assignments, effective n = number of clusters, CV folds at the cluster level).

#' Combined randomization-based test of the sharp null
#'
#' @param formula a design formula `Y ~ treatment` for complete randomization,
#'   or `Y ~ treatment | block` for randomization within blocks (or use `blocks`).
#' @param data a data.frame/data.table.
#' @param blocks optional block specification: a vector of labels or a column
#'   name in `data`. Omit when the block is in the formula as `| block`.
#'   With neither specification, the design fixes the total treated count.
#' @param clusters optional cluster specification for a cluster-randomized
#'   design: a vector of cluster labels or a column name in `data`. When given,
#'   the design is collapsed to one row per cluster (outcome and covariates
#'   aggregated by `cluster_agg`), treatment is permuted over clusters within
#'   blocks, and the effective sample size is the number of clusters. Treatment
#'   must be constant within a cluster and clusters must nest within blocks.
#' @param statistic which combination to use: `"screen"` (default), `"quadratic"`,
#'   `"cauchy"`, or `"max"`.
#' @param representations a named list of representation functions;
#'   [riposte_reps_default()] by default. Supply your own to combine a different
#'   set.
#' @param nresample number of re-randomizations for the Monte-Carlo combinations.
#'   Ignored for `engine = "asymptotic"`.
#' @param screen a [riposte_screen_control()] list, used when
#'   `statistic = "screen"`.
#' @param adjust optional covariance adjustment: either a one-sided covariate
#'   formula (e.g. `~ x1 + x2`) or a [riposte_adjust()] specification. When a bare
#'   formula, it is wrapped with `learner` (ridge by default). `NULL` means no
#'   adjustment. The combination is then computed on the controls-only residuals,
#'   refit inside each re-randomization; see [riposte_adjust()] for what is and is
#'   not exact under adjustment.
#' @param learner optional controls-only learner `(X, Y, control) -> Yhat` for
#'   covariance adjustment; [riposte_ridge_learner()] by default. Overrides the
#'   learner in an `adjust` specification if both are given.
#' @param adjust_refit refit the learner inside every re-randomization (`TRUE`,
#'   exact for any learner) or fit once and permute fixed residuals (`FALSE`, the
#'   shortcut). Used only when `adjust` is a bare formula.
#' @param cluster_agg how to aggregate the outcome (and covariates) within a
#'   cluster when `clusters` is given; the cluster mean by default.
#' @param seed optional integer seed (L'Ecuyer-CMRG) for reproducibility.
#'   Ignored for `engine = "asymptotic"`, which uses no random draws.
#' @param engine `"permute"` (default) refers every combination to brute-force
#'   re-randomization (exact). `"asymptotic"` uses a chi-square reference for
#'   `statistic = "quadratic"`, with degrees of freedom equal to the numerical
#'   rank of the permutation covariance. For `statistic = "cauchy"`, it computes
#'   each representation's two-sided p-value from its squared standardized
#'   statistic and a chi-square distribution on one degree of freedom, then
#'   compares the mean Cauchy transform with a standard Cauchy distribution.
#'   These are approximations, available without `adjust`, and do not require
#'   `coin` or `fastperm`. Set `statistic` explicitly to `"quadratic"` or
#'   `"cauchy"`; the max and screen approximations are not implemented.
#'   `"saddlepoint"` computes the combination with no
#'   re-randomization, via `fastperm`: the unadjusted Cauchy (each representation's
#'   marginal saddlepoint p-value combined by the Liu-Xie analytic tail), the
#'   unadjusted quadratic (the omnibus `Q` referred to fastperm's multivariate
#'   saddlepoint M2, feasible for a few representations), or `"screen"`, which uses
#'   the permutation `Sigma` to choose between them (and to regularize a
#'   rank-deficient metric). It is an approximation, supports `statistic` other than
#'   `"max"` without `adjust`, and requires the `fastperm` package (the quadratic
#'   needs a `fastperm` with `fastperm_spa_quadratic`).
#' @param ... reserved.
#' @return an object with the observed statistic, the p-value, the
#'   combination used, the condition number, and (for `"screen"`) the choice made.
#'   With `engine = "permute"`, the p-value is the share of the
#'   `nresample` re-randomizations at least as extreme as observed, with the
#'   observed assignment included. That observed-inclusion convention makes the
#'   test's level exactly valid in finite samples for any `nresample` (the p-value
#'   estimates, rather than equals, the full-enumeration value, so it varies with
#'   the seed). The closed-form moments and the screen's condition number are
#'   exact, not Monte-Carlo. With `engine = "asymptotic"`, `nresample` is zero
#'   and the p-value uses the stated approximation, not a permutation count.
#'   Cauchy results from this engine also include `component_p`, the individual
#'   chi-square p-values.
#'
#' @details
#' The chi-square approximation requires a normal limit for the standardized
#' treated-score sums on their nondegenerate subspace. A large number of units
#' alone does not ensure this when individual scores dominate their variance.
#' For the Cauchy combination, dependent component p-values do not in general
#' yield an exactly standard Cauchy statistic: the analytic upper tail is an
#' approximation, not a finite-sample level guarantee.
#'
#' The asymptotic Cauchy calculation retains both tails on the log scale to
#' avoid artificial zeros and ones. An exactly zero component statistic gives
#' an individual p-value of one and a Cauchy term of negative infinity; the
#' combined p-value is then one unless an opposing infinite term makes the sum
#' undefined. That case gives an error. This engine never switches silently to
#' permutations. It combines the individual p-values only; it does not add
#' the quadratic p-value as another term. The example in [riposte_acat_term()]
#' shows how to include that p-value to calculate the paper's Hybrid test.
#' @export
riposte_test <- function(formula, data, blocks = NULL, clusters = NULL,
                        statistic = c("screen", "quadratic", "cauchy", "max"),
                        representations = riposte_reps_default(),
                        nresample = 1999L,
                        screen = riposte_screen_control(),
                        adjust = NULL, learner = NULL, adjust_refit = TRUE,
                        cluster_agg = mean, seed = NULL,
                        engine = c("permute", "saddlepoint", "asymptotic"), ...) {
  statistic <- match.arg(statistic)
  engine <- match.arg(engine)
  if (engine == "asymptotic") {
    if (!statistic %in% c("quadratic", "cauchy"))
      stop("engine = \"asymptotic\" requires statistic = \"quadratic\" or ",
           "\"cauchy\". Use engine = \"permute\" for the max or screen.", call. = FALSE)
    if (!is.null(adjust))
      stop("engine = \"asymptotic\" does not support covariance adjustment; ",
           "use engine = \"permute\".", call. = FALSE)
    nresample <- 0L
  } else {
    nresample <- riposte_check_nresample(nresample)
    if (!is.null(seed)) { RNGkind("L'Ecuyer-CMRG"); set.seed(seed) }
  }

  ## the saddlepoint engine refers the Cauchy, the quadratic, and the screen's
  ## choice between them to fastperm without re-randomizing; the max needs a joint
  ## saddlepoint not yet built, and adjustment is not available draws-free
  if (engine == "saddlepoint") {
    if (statistic == "max")
      stop("engine = \"saddlepoint\" does not support statistic = \"max\"; the max ",
           "needs a joint saddlepoint that is not yet implemented. Use ",
           "engine = \"permute\".", call. = FALSE)
    if (!is.null(adjust))
      stop("engine = \"saddlepoint\" does not support covariance adjustment; ",
           "use engine = \"permute\".", call. = FALSE)
  }

  des <- riposte_parse_design(formula, data, blocks)

  ## build the (unit-level) covariate matrix once if adjusting, so it can be
  ## collapsed alongside the design under cluster randomization
  spec <- NULL; X <- NULL
  if (!is.null(adjust)) {
    spec <- if (inherits(adjust, "riposte_adjust")) adjust
            else riposte_adjust(adjust,
                                learner = learner %||% riposte_ridge_learner(),
                                refit = adjust_refit)
    if (!is.null(learner)) spec$learner <- learner
    X <- riposte_covariate_matrix(spec$covariates, data)
  }

  ## cluster randomization: collapse to one row per cluster (the sampling unit)
  n_units <- length(des$y); n_clusters <- NA_integer_
  if (!is.null(clusters)) {
    cl <- riposte_resolve_clusters(clusters, data, n_units)
    col <- riposte_collapse_clusters(des$y, des$z, des$block, cl, X, cluster_agg)
    des <- list(y = col$y, z = col$z, block = col$block)
    X <- col$X; n_clusters <- col$n_clusters
  }

  ## Only the permutation engine needs randomized assignments.
  draws <- if (engine == "permute")
    riposte_block_draws(des$z, des$block, nresample) else NULL

  out <- if (is.null(adjust))
    riposte_test_unadjusted(des, representations, draws, statistic, screen, engine)
  else
    riposte_test_adjusted(des, X, representations, draws, statistic, screen, spec)

  structure(c(out, list(
    requested = statistic,
    engine = engine,
    n = length(des$y),
    n_units = n_units,
    n_clusters = n_clusters,
    nblocks = nlevels(des$block),
    nresample = nresample,
    adjusted = !is.null(adjust),
    clustered = !is.null(clusters),
    call = match.call()
  )), class = "riposte_test")
}

## unadjusted path: closed-form moments + the shared draws
riposte_test_unadjusted <- function(des, representations, draws, statistic, screen,
                                    engine = "permute") {
  sm <- riposte_score_matrix(des$y, des$block, representations)
  if (ncol(sm$scores) == 0L)
    stop("every representation is constant within blocks; nothing to test.",
         call. = FALSE)

  if (engine == "asymptotic")
    return(riposte_test_asymptotic(des, sm, statistic))

  ## draws-free saddlepoint combinations. The Cauchy needs no joint covariance; the
  ## quadratic (and the screen's choice between it and the Cauchy) need the
  ## permutation Sigma, which also lets the screen regularize a rank-deficient metric
  ## (shrink) or fall back to the Cauchy (threshold) instead of inverting a singular M.
  if (engine == "saddlepoint") {
    chosen <- statistic
    scr <- NULL
    metric <- "cov"
    if (statistic %in% c("screen", "quadratic")) {
      moments <- riposte_sw_moments(sm$scores, des$z, des$block)
      riposte_assert_testable(moments$Sigma)
      scr <- riposte_screen(moments$Sigma, screen)
      if (statistic == "screen") {
        chosen <- if (scr$method == "threshold") scr$choice else "quadratic"
        if (scr$method == "shrink") metric <- scr$Sigma
      }
    }
    res <- if (chosen == "quadratic")
      riposte_quadratic_spa(sm$scores, des$z, des$block, metric = metric)
    else
      riposte_cauchy_spa(sm$scores, des$z, des$block)
    return(list(
      statistic = res$statistic, p.value = res$p.value, combination = chosen,
      condition = if (is.null(scr)) NA_real_ else scr$condition,
      screen = if (statistic == "screen") scr else NULL,
      df = if (chosen == "quadratic") res$rank else NA_integer_,
      kept = sm$kept, dropped = sm$dropped))
  }

  moments <- riposte_sw_moments(sm$scores, des$z, des$block)
  riposte_assert_testable(moments$Sigma)
  scr <- riposte_screen(moments$Sigma, screen)
  chosen <- statistic
  if (statistic == "screen")
    chosen <- if (scr$method == "threshold") scr$choice else "quadratic"

  res <- switch(
    chosen,
    quadratic = {
      mom <- if (statistic == "screen" && scr$method == "shrink")
        list(mu = moments$mu, Sigma = scr$Sigma) else moments
      riposte_quadratic(sm$scores, des$z, des$block, moments = mom, draws = draws)
    },
    cauchy = riposte_cauchy(sm$scores, des$z, des$block, draws = draws),
    max    = riposte_max(sm$scores, des$z, des$block, moments = moments, draws = draws)
  )
  list(statistic = res$statistic, p.value = res$p.value, combination = chosen,
       condition = scr$condition, screen = if (statistic == "screen") scr else NULL,
       df = res$df, kept = sm$kept, dropped = sm$dropped)
}

## adjusted path: refit-per-permutation residuals + Monte-Carlo moments. X is the
## (possibly cluster-collapsed) covariate matrix built by the caller.
riposte_test_adjusted <- function(des, X, representations, draws, statistic,
                                  screen, spec) {
  rp <- riposte_residual_perm_stats(des$y, X, des$block, draws, spec$learner,
                                    representations, spec$refit)
  comb <- riposte_adjusted_combination(rp$stats, statistic, screen)
  list(statistic = comb$result$statistic, p.value = comb$result$p.value,
       combination = comb$chosen, condition = comb$condition,
       screen = comb$screen, df = comb$result$df, kept = rp$kept,
       dropped = setdiff(names(representations), rp$kept),
       refit = spec$refit)
}
