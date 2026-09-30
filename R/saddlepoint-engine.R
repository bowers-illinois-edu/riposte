## saddlepoint-engine.R --- the optional draws-free Cauchy via fastperm.
##
## THE POINT. The default engine refers every combination to brute-force
## re-randomization. For the UNADJUSTED Cauchy combination there is a closed-form
## alternative: each representation's two-sided permutation p-value comes from a
## saddlepoint approximation to the exact within-block permutation distribution
## (no draws), and the Liu-Xie analytic tail combines them. We verified this
## controls the false-positive rate (dev/spike-acat-size*.R) and that the
## saddlepoint marginals match the brute-force mid-p marginals to ~1e-3
## (dev/spike-acat-fast.R). The saddlepoint lives in the separate `fastperm`
## package; riposte reaches it via requireNamespace (an OPTIONAL runtime
## dependency, not Suggests/Remotes), so a missing or non-CRAN fastperm cannot
## break riposte's own checks or CI.
##
## This is an APPROXIMATION engine. The default permutation Cauchy stays exact;
## this trades the exact permutation calibration for a closed-form, draws-free
## analytic combination whose p-value differs from the permutation one by the
## analytic-vs-permutation gap (about 0.02 in the body, tighter in the tail).

#' Draws-free Cauchy combination via the fastperm saddlepoint
#'
#' Computes the unadjusted Cauchy (ACAT) combination without re-randomizing: each
#' representation's two-sided within-block permutation p-value is the fastperm
#' saddlepoint p-value, and the Liu-Xie analytic standard-Cauchy tail combines
#' them. Requires the `fastperm` package.
#'
#' @param scores within-block-centred score matrix (units x representations).
#' @param z 0/1 treatment vector.
#' @param block block factor.
#' @return a list with the ACAT `statistic`, its analytic `p.value`, and
#'   `component_p` (each representation's two-sided saddlepoint p-value).
#' @keywords internal
#' @noRd
riposte_cauchy_spa <- function(scores, z, block) {
  ## fastperm is an OPTIONAL dependency (Suggests + Remotes: it lives on GitHub,
  ## not CRAN); guard so riposte works without it and only this engine needs it
  if (!requireNamespace("fastperm", quietly = TRUE))
    stop("engine = \"saddlepoint\" requires the fastperm package. Install it with\n",
         "  remotes::install_github(\"bowers-illinois-edu/fastperm\")", call. = FALSE)
  scores <- as.matrix(scores)

  ## two-sided saddlepoint permutation p-value of each representation's linear
  ## statistic, computed from the exact permutation CGF with no draws
  p <- vapply(seq_len(ncol(scores)), function(j)
    fastperm::fastperm_spa_linear(scores[, j], z, block,
                                  alternative = "two.sided")$p.value, numeric(1))
  names(p) <- colnames(scores)

  ## Liu-Xie analytic Cauchy combination: the ACAT statistic is the mean of the
  ## pole-aware Cauchy transforms, and its tail is approximated by a standard
  ## Cauchy, so the combined p-value is 0.5 - atan(T)/pi -- no permutation.
  Tc <- mean(riposte_acat_term(p))
  list(statistic = unname(Tc), p.value = 0.5 - atan(Tc) / pi, component_p = p)
}

#' Draws-free quadratic (omnibus) combination via the fastperm saddlepoint
#'
#' Computes the Hansen-Bowers omnibus quadratic
#' `Q = (T - mu)' Sigma^{-1} (T - mu)` without re-randomizing: `T` is the vector of
#' within-block treated-score sums, `mu` and `Sigma` are its exact permutation mean
#' and covariance, and the within-block permutation tail of `Q` is fastperm's
#' multivariate saddlepoint (M2), which carries the true skewness and higher
#' cumulants of the permutation law rather than a Gaussian approximation. Requires a
#' fastperm that has `fastperm_spa_quadratic` (the route-b feature).
#'
#' @param scores within-block-centred score matrix (units x representations). `Q` is
#'   invariant to within-block centring, so this is the same `Q` as the permutation
#'   quadratic's.
#' @param z 0/1 treatment vector.
#' @param block block factor.
#' @param metric `"cov"` (the exact permutation covariance of `T`, the Hansen-Bowers
#'   omnibus) or a symmetric positive-definite matrix. The screen passes its shrunk,
#'   well-conditioned covariance here to regularize a rank-deficient metric.
#' @return a list with the omnibus `statistic` (`Q`), its saddlepoint `p.value`
#'   (Lugannani-Rice), and the metric `rank` (the chi-square degrees of freedom).
#' @keywords internal
#' @noRd
riposte_quadratic_spa <- function(scores, z, block, metric = "cov") {
  if (!requireNamespace("fastperm", quietly = TRUE))
    stop("engine = \"saddlepoint\" requires the fastperm package. Install it with\n",
         "  remotes::install_github(\"bowers-illinois-edu/fastperm\")", call. = FALSE)
  ## the quadratic path needs a newer fastperm than the scalar Cauchy path does
  if (!exists("fastperm_spa_quadratic", where = asNamespace("fastperm")))
    stop("statistic = \"quadratic\" with engine = \"saddlepoint\" needs a fastperm ",
         "that has fastperm_spa_quadratic (the multivariate saddlepoint); the ",
         "installed version lacks it. Update fastperm, or use engine = \"permute\".",
         call. = FALSE)
  scores <- as.matrix(scores)

  ## metric = "cov" makes fastperm use the exact permutation covariance of T, so Q
  ## is the Hansen-Bowers omnibus -- identical to the permutation quadratic's Q (a
  ## matrix metric, e.g. the screen's shrunk Sigma, regularizes a rank-deficient
  ## set). M2 inverts the exact permutation CGF of Q; its tensor Gauss-Hermite grid
  ## is feasible only for a few representations, so a high-rank set errors -- re-raise
  ## that one with the riposte-level way out, but let other errors pass through.
  res <- tryCatch(
    fastperm::fastperm_spa_quadratic(scores, z, block, metric = metric,
                                     method = "saddlepoint"),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("tensor|nodes|sparse-grid|representations", msg))
        stop("the quadratic saddlepoint uses a tensor grid feasible only for a few ",
             "representations (", msg, "). Use engine = \"permute\" (exact for any ",
             "number of representations), or supply fewer representations.",
             call. = FALSE)
      stop(e)
    })
  list(statistic = res$statistic, p.value = res$p.value, rank = res$rank)
}
