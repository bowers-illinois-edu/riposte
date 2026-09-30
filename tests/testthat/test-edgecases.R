## test-edgecases.R
##
## Regression tests for the degenerate-input handling found in the correctness
## audit. Each asserts a CLEAR, specific error (or a clean result) instead of the
## opaque downstream failure the audit reproduced.

make_df <- function(seed = 1) {
  set.seed(seed)
  block <- factor(rep(1:4, each = 6)); N <- 24
  z <- integer(N)
  for (b in levels(block)) { ix <- which(block == b); z[ix][sample.int(6, 3)] <- 1L }
  data.frame(Y = rnorm(N), trt = z, blk = block)
}

test_that("a single-observation score matrix does not crash", {
  ## n == 1: every representation is constant within the lone block, so the
  ## result is a clean 1-by-0 matrix, not an opaque colnames error
  sm <- riposte_score_matrix(c(5), factor(1))
  expect_equal(nrow(sm$scores), 1L)
  expect_equal(ncol(sm$scores), 0L)
})

test_that("missing values in treatment, outcome, or block are rejected clearly", {
  d <- make_df()
  d_trt <- d; d_trt$trt[1] <- NA
  expect_error(riposte_test(Y ~ trt | blk, d_trt, nresample = 99), "treatment has 1 missing")
  d_y <- d; d_y$Y[1] <- NA
  expect_error(riposte_test(Y ~ trt | blk, d_y, nresample = 99), "outcome has 1 missing")
  d_b <- d; d_b$blk[1] <- NA
  expect_error(riposte_test(Y ~ trt | blk, d_b, nresample = 99), "block has 1 missing")
})

test_that("nresample must be a single whole number >= 1", {
  d <- make_df()
  expect_error(riposte_test(Y ~ trt | blk, d, nresample = -5), "whole number")
  expect_error(riposte_test(Y ~ trt | blk, d, nresample = 10.7), "whole number")
  expect_error(riposte_test(Y ~ trt | blk, d, nresample = 0), "whole number")
  expect_error(riposte_components(Y ~ trt | blk, d, nresample = 0), "whole number")
})

test_that("a design with no within-block randomization errors with a named cause", {
  ## every block all-control: the Strasser-Weber weight is zero everywhere, so
  ## there is nothing to test -- and the user is told exactly that
  dd <- data.frame(Y = rnorm(8), trt = rep(0, 8), blk = rep(c("a", "b"), each = 4))
  expect_error(riposte_test(Y ~ trt | blk, dd, nresample = 99),
               "no within-block randomization")
  da <- data.frame(Y = rnorm(8), trt = rep(1, 8), blk = rep(c("a", "b"), each = 4))
  expect_error(riposte_test(Y ~ trt | blk, da, nresample = 99),
               "no within-block randomization")
})

test_that("a partially degenerate design (one block degenerate) still works", {
  ## only one block all-control: the other block carries the test
  d <- data.frame(Y = rnorm(12), trt = c(rep(0, 6), rep(c(1, 0), 3)),
                  blk = rep(c("a", "b"), each = 6))
  r <- riposte_test(Y ~ trt | blk, d, statistic = "max", nresample = 99, seed = 1)
  expect_true(r$p.value > 0 && r$p.value <= 1)
})
