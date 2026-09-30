## learners.R --- controls-only learners for covariance adjustment.
##
## THE POINT. A learner predicts the outcome from covariates using the CONTROL
## units only, so the prediction never sees the treated outcomes and cannot
## absorb the treatment effect. riposte residualizes Y by this prediction and
## tests the residuals. The contract is deliberately minimal so any method fits:
##
##   learner(X, Y, control)  ->  Yhat   (a length-N vector of fitted values)
##
## where `control` is a logical vector marking the units to fit on. The same
## learner is called afresh inside every re-randomization (refit-per-permutation),
## with `control` set to that draw's control units; this is what keeps the test
## exact for any learner (the statistic is then a deterministic function of the
## assignment, the fixed outcomes, and the covariates). A learner that is cheap to
## refit (ridge has a closed form) keeps the permutation loop fast.

#' Controls-only ridge learner
#'
#' Returns a learner that fits ridge regression on the control units and predicts
#' for all units. Ridge has a closed form (no iterative fitting) and its
#' shrinkage keeps predictions stable when there are many covariates --- the
#' regime where an unregularized fit overfits and the fit-once shortcut loses
#' exactness.
#'
#' @param lambda ridge penalty. If `NULL`, chosen by leave-one-out CV on the
#'   control units each time the learner is called. A fixed `lambda` (for
#'   instance the CV value from the observed controls) is faster and equally
#'   exact under refit-per-permutation.
#' @return a function `(X, Y, control) -> Yhat`.
#' @export
riposte_ridge_learner <- function(lambda = NULL) {
  function(X, Y, control) {
    X <- as.matrix(X)
    Xc <- X[control, , drop = FALSE]
    Yc <- Y[control]
    xm <- colMeans(Xc); ym <- mean(Yc)
    Xcc <- sweep(Xc, 2, xm)               # centre so the intercept is unpenalized
    Ycc <- Yc - ym
    lam <- if (is.null(lambda)) riposte_ridge_cv(Xcc, Ycc) else lambda
    G <- crossprod(Xcc) + lam * diag(ncol(X))
    beta <- solve(G, crossprod(Xcc, Ycc))
    as.numeric(sweep(X, 2, xm) %*% beta) + ym
  }
}

#' Controls-only ordinary-least-squares learner
#'
#' Returns a learner that fits OLS on the control units. With many covariates
#' relative to the number of controls, OLS overfits --- which is exactly the case
#' that makes the fit-once / fixed-residual shortcut anti-conservative and that
#' refit-per-permutation handles correctly. Useful for demonstrating that
#' contrast; ridge is the better default in practice.
#'
#' @return a function `(X, Y, control) -> Yhat`.
#' @export
riposte_lm_learner <- function() {
  function(X, Y, control) {
    X <- as.matrix(X)
    Xc <- X[control, , drop = FALSE]
    Yc <- Y[control]
    xm <- colMeans(Xc); ym <- mean(Yc)
    Xcc <- sweep(Xc, 2, xm)
    Ycc <- Yc - ym
    ## minimum-norm OLS via the pseudoinverse (handles d >= n_control)
    beta <- MASS::ginv(crossprod(Xcc)) %*% crossprod(Xcc, Ycc)
    as.numeric(sweep(X, 2, xm) %*% beta) + ym
  }
}

#' Leave-one-out CV ridge penalty via the SVD shortcut
#'
#' Picks the ridge penalty minimizing the leave-one-out CV error, evaluated for
#' all candidate `lambda` at once from one SVD. The learner fits ridge with an
#' UNPENALIZED intercept (it centres on the control means), so the hat-matrix
#' diagonal carries the intercept term `1 / n_control` on top of the penalized
#' part: `h_ii = 1/n_control + sum_j u_{ij}^2 s_j^2 / (s_j^2 + lambda)`, and
#' `LOOCV(lambda) = mean( ((Yc - Yhat_i) / (1 - h_ii))^2 )`. Dropping the
#' `1/n_control` term understates leverage and makes CV pick a too-small penalty
#' exactly in the high-dimensional regime ridge exists for (where it would
#' interpolate the controls).
#'
#' @param Xcc centred control design matrix (n_control x d).
#' @param Ycc centred control outcome.
#' @param lambda_seq candidate penalties; a log-spaced grid by default.
#' @return the minimizing `lambda` (a positive scalar).
#' @keywords internal
#' @noRd
riposte_ridge_cv <- function(Xcc, Ycc,
                             lambda_seq = 10^seq(-3, 3, length.out = 50)) {
  sv <- svd(Xcc)
  U <- sv$u; s <- sv$d
  UtY <- crossprod(U, Ycc)                # U' Yc
  U2 <- U^2
  n_control <- length(Ycc)
  errs <- vapply(lambda_seq, function(lam) {
    filt <- s^2 / (s^2 + lam)
    Yhat <- as.numeric(U %*% (filt * UtY))
    ## leverage h_ii: the unpenalized intercept contributes 1/n_control on top of
    ## the penalized SVD part; omitting it understates leverage and biases CV
    ## toward too-small penalties (verified against brute-force LOOCV).
    h <- 1 / n_control + as.numeric(U2 %*% filt)
    r <- (Ycc - Yhat) / (1 - h)
    mean(r^2)
  }, numeric(1))
  lambda_seq[which.min(errs)]
}
