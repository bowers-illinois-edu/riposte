## score.R --- turn representations into within-block-centred score columns.
##
## THE POINT. Each representation becomes a Strasser-Weber linear statistic: at
## an assignment z, the statistic for representation c is sum_i z_i * s_{c,i},
## where s_{c,i} is unit i's score CENTRED within its block. Centring within
## block makes each statistic have permutation mean zero (the treated count per
## block is fixed), which is what the closed-form moments and the combinations
## assume. This function applies the representation functions block by block and
## returns the centred score matrix, dropping any representation that is constant
## in every block (it carries no information and would be a zero column).

#' Build the within-block-centred score matrix from representations
#'
#' Applies each representation function to the outcome within each block, then
#' centres each resulting score within block (subtract the block mean). Returns
#' the matrix of centred scores, one column per surviving representation. A
#' representation that is constant within every block contributes nothing to a
#' permutation statistic and is dropped, with a recorded reason.
#'
#' @param y numeric outcome vector.
#' @param block factor or vector of block labels, the same length as `y`.
#' @param representations a named list of representation functions; see
#'   [riposte_reps_default()].
#' @return a list with
#'   \describe{
#'     \item{scores}{numeric matrix, `length(y)` rows, one column per kept
#'       representation, centred within block.}
#'     \item{kept}{character vector of representation names retained.}
#'     \item{dropped}{character vector of representation names dropped for being
#'       within-block constant.}
#'   }
#' @export
riposte_score_matrix <- function(y, block, representations = riposte_reps_default()) {
  centred <- riposte_centred_scores(y, block, representations)
  ## drop representations that are constant within every block (all-zero column)
  keep <- colSums(abs(centred)) > 1e-12
  nm <- colnames(centred)
  list(
    scores  = centred[, keep, drop = FALSE],
    kept    = nm[keep],
    dropped = nm[!keep]
  )
}

#' Within-block-centred scores for every representation (no dropping)
#'
#' The primitive behind [riposte_score_matrix()]: applies each representation
#' within block and centres it, returning one column per representation in the
#' order given, WITHOUT dropping constant columns. The covariance-adjustment
#' engine needs this no-drop form so the column set stays the same across
#' re-randomizations (a representation that is constant within a block on one
#' draw's residuals contributes a zero column, not a missing one).
#'
#' @inheritParams riposte_score_matrix
#' @return a numeric matrix, `length(y)` rows, one column per representation.
#' @keywords internal
#' @noRd
riposte_centred_scores <- function(y, block, representations = riposte_reps_default()) {
  riposte_validate_reps(representations)
  y <- as.numeric(y)
  n <- length(y)
  if (length(block) != n)
    stop("`y` and `block` must have the same length.", call. = FALSE)
  block <- as.factor(block)
  idx_list <- split(seq_len(n), block)

  ## apply one representation across all blocks, returning a length-n vector in
  ## the original row order; check the function returns the right length per block
  apply_rep <- function(f, nm) {
    out <- numeric(n)
    for (ix in idx_list) {
      val <- f(y[ix])
      if (length(val) != length(ix))
        stop(sprintf(
          "representation '%s' returned length %d for a block of size %d; representations must return one score per unit.",
          nm, length(val), length(ix)), call. = FALSE)
      out[ix] <- as.numeric(val)
    }
    out
  }

  nm <- names(representations)
  ## matrix() keeps the k columns even when n == 1, where vapply would otherwise
  ## simplify to a length-k vector and the colnames assignment would fail
  raw <- matrix(
    vapply(seq_along(representations),
           function(j) apply_rep(representations[[j]], nm[j]), numeric(n)),
    nrow = n, dimnames = list(NULL, nm))

  ## centre each score within block: subtract the within-block mean so the
  ## permutation mean of the linear statistic is zero
  centred <- apply(raw, 2, function(col) col - stats::ave(col, block))
  if (is.null(dim(centred)))
    centred <- matrix(centred, nrow = n, dimnames = list(NULL, nm))
  centred
}
