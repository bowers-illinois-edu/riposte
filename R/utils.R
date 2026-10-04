## utils.R --- parse the user's design into outcome, treatment, and block.
##
## The user describes the design with a formula. We accept the following
## spellings and reduce them to three vectors the rest of the package
## uses: the outcome y, a 0/1 treatment z, and a block factor.
##   - coin style:   Y ~ treatment | block        (block in the formula)
##   - blocks arg:   Y ~ treatment, blocks = ...   (block given separately)
##   - complete:     Y ~ treatment                 (one randomized group)
## Giving the block both ways is ambiguous and is an error.

#' Parse a design formula into outcome, treatment, and block vectors
#'
#' @param formula a formula `Y ~ treatment` or `Y ~ treatment | block`.
#' @param data a data.frame or data.table holding the variables.
#' @param blocks optional block specification when the formula has no `| block`:
#'   either a vector of block labels (length matching `data`) or a single column
#'   name in `data`. With neither specification, use complete randomization:
#'   one group with the observed total treated count fixed.
#' @return a list with `y` (numeric), `z` (0/1 numeric treatment), and `block`
#'   (factor).
#' @keywords internal
#' @noRd
riposte_parse_design <- function(formula, data, blocks = NULL) {
  if (!inherits(formula, "formula"))
    stop("`formula` must be a formula, e.g. Y ~ treatment | block.", call. = FALSE)
  data <- as.data.frame(data)
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]

  ## a `| block` term parses as a call to `|`; split it out if present
  block_expr <- NULL
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    trt_expr   <- rhs[[2L]]
    block_expr <- rhs[[3L]]
  } else {
    trt_expr <- rhs
  }

  if (!is.null(block_expr) && !is.null(blocks))
    stop("give the block either in the formula (`| block`) or via `blocks`, not both.",
         call. = FALSE)

  y <- as.numeric(eval(lhs, data, environment(formula)))
  z <- riposte_as_indicator(eval(trt_expr, data, environment(formula)))

  block <- if (!is.null(block_expr)) {
    eval(block_expr, data, environment(formula))
  } else if (!is.null(blocks)) {
    if (length(blocks) == 1L && is.character(blocks) && blocks %in% names(data))
      data[[blocks]]
    else blocks
  } else {
    rep(1L, length(y))
  }

  if (length(block) != length(y))
    stop("block length does not match the data.", call. = FALSE)

  ## a randomization test conditions on the realized design; silently dropping
  ## units with missing values would change it, so reject them with a clear,
  ## variable-specific message rather than failing later in the linear algebra
  for (v in list(c("outcome", "y"), c("treatment", "z"), c("block", "block"))) {
    val <- get(v[2])
    if (anyNA(val))
      stop(sprintf("%s has %d missing value(s); remove or impute them before testing.",
                   v[1], sum(is.na(val))), call. = FALSE)
  }

  list(y = y, z = z, block = as.factor(block))
}

#' Validate the re-randomization count
#'
#' `nresample` must be a single whole number at least 1. A negative value would
#' silently collapse the p-value to 1, and a non-integer would break the printed
#' summary; both are rejected up front.
#'
#' @param nresample the value to check.
#' @return `nresample` as an integer, invisibly, if valid.
#' @keywords internal
#' @noRd
riposte_check_nresample <- function(nresample) {
  if (length(nresample) != 1L || !is.finite(nresample) ||
      nresample < 1 || nresample %% 1 != 0)
    stop("`nresample` must be a single whole number >= 1.", call. = FALSE)
  as.integer(nresample)
}

#' Stop if no block has any within-block randomization
#'
#' When every block is all-treated or all-control, the Strasser-Weber weight is
#' zero in every block, so the permutation covariance is all zeros and the
#' downstream `cov2cor`/`eigen` fail with an opaque message. Catch that here and
#' name the real cause.
#'
#' @param Sigma a permutation covariance matrix.
#' @keywords internal
#' @noRd
riposte_assert_testable <- function(Sigma) {
  ## a design with no block holding both arms gives every variance exactly
  ## zero (each block's weight m (n - m) / (n (n - 1)) is 0), so the check is
  ## for exact zeros; a fixed cutoff such as 1e-12 refused testable designs
  ## when Y was recorded in small units (issue #1)
  if (length(Sigma) == 0L || all(diag(Sigma) <= 0))
    stop("no block has both treated and control units; there is no within-block ",
         "randomization to test.", call. = FALSE)
}

#' Coerce a treatment variable to a 0/1 indicator
#'
#' Numeric 0/1 stays as is; a two-level factor/character/logical maps its second
#' level (or TRUE) to 1. Anything else is an error: riposte tests a binary
#' treatment.
#'
#' @param x a treatment vector.
#' @return a 0/1 numeric vector.
#' @keywords internal
#' @noRd
riposte_as_indicator <- function(x) {
  if (is.logical(x)) return(as.numeric(x))
  if (is.numeric(x)) {
    u <- unique(x[!is.na(x)])
    if (all(u %in% c(0, 1))) return(as.numeric(x))
    stop("numeric treatment must be coded 0/1.", call. = FALSE)
  }
  f <- as.factor(x)
  if (nlevels(f) != 2L)
    stop("treatment must have exactly two levels.", call. = FALSE)
  as.numeric(f == levels(f)[2L])
}

#' Moore-Penrose inverse of a covariance matrix on the correlation scale
#'
#' The quadratic Q = c' Sigma^+ c is unchanged, in exact arithmetic, when each
#' score sum is divided by its standard deviation (c' Sigma^+ c =
#' (D^-1 c)' R^+ (D^-1 c) for c in the column space of Sigma, with
#' R = D^-1 Sigma D^-1). In floating point the two differ: the eigenvalues of
#' Sigma carry the units of the scores and can span ten or more orders of
#' magnitude, so a cutoff relative to the largest drops real directions, while
#' the eigenvalues of R sum to the number of scores whatever the units. The
#' cutoff is MASS::ginv()'s, sqrt(eps) times the largest eigenvalue, and the
#' rank counts the directions kept, so the degrees of freedom match Q.
#'
#' @param Sigma a covariance matrix of score sums.
#' @return a list with `sd` (the standard deviations, with 1 in place of 0 so
#'   that a zero-variance sum contributes nothing), `R_pinv` (the Moore-Penrose
#'   inverse of the correlation matrix), and `rank` (the directions kept).
#' @keywords internal
#' @noRd
riposte_std_pinv <- function(Sigma) {
  sdv <- sqrt(pmax(diag(Sigma), 0))
  sdv[sdv == 0] <- 1
  R <- Sigma / outer(sdv, sdv)
  e <- eigen(R, symmetric = TRUE)
  keep <- e$values > max(e$values, 0) * sqrt(.Machine$double.eps)
  V <- e$vectors[, keep, drop = FALSE]
  R_pinv <- V %*% (t(V) / e$values[keep])
  dimnames(R_pinv) <- dimnames(Sigma)
  list(sd = sdv, R_pinv = R_pinv, rank = sum(keep))
}
