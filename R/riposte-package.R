#' riposte: combined randomization-based tests across outcome representations
#'
#' riposte provides tests of the sharp null that combine several
#' representations of an outcome to detect treatment effects a difference of
#' means misses, in completely randomized, block-randomized, and
#' cluster-randomized experiments. By default, p-values use direct
#' re-randomization under the specified design. For unadjusted quadratic and
#' Cauchy tests, users can instead choose `engine = "asymptotic"` to calculate
#' p-values with chi-square and Cauchy approximations. See [riposte_test()].
#'
#' The three combinations are:
#' \itemize{
#'   \item the quadratic (energy/distance) omnibus, which uses the full
#'         permutation covariance of the representations;
#'   \item the Cauchy combination, which uses only each representation's
#'         marginal calibration; and
#'   \item the max combination.
#' }
#' With the permutation engine, a screen chooses between the quadratic and
#' the Cauchy from the conditioning
#' of the permutation covariance, an ancillary quantity (a function of the
#' outcomes under the sharp null and the design, not of the realized
#' assignment).
#'
#' @keywords internal
"_PACKAGE"
