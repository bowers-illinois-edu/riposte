## test-integration.R
##
## End-to-end behaviour of the user-facing functions: the design parses, the
## chosen combination runs, the result is reproducible under a seed, the screen
## reports its choice, and the not-yet-built paths error honestly.

make_df <- function(seed = 1, effect = 0) {
  set.seed(seed)
  block <- factor(rep(1:6, each = 10))
  z <- integer(length(block))
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(10, 5)] <- 1L
  }
  y <- rnorm(length(block)) + effect * z
  data.frame(Y = y, trt = z, blk = block)
}

test_that("riposte_test runs each combination and returns a valid p-value", {
  d <- make_df()
  for (st in c("quadratic", "cauchy", "max", "screen")) {
    r <- riposte_test(Y ~ trt | blk, d, statistic = st, nresample = 199, seed = 11)
    expect_s3_class(r, "riposte_test")
    expect_true(r$p.value > 0 && r$p.value <= 1)
    expect_output(print(r), "riposte test")
  }
})

test_that("the formula `| block` and the `blocks` argument agree", {
  d <- make_df()
  a <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic", nresample = 199, seed = 5)
  b <- riposte_test(Y ~ trt, d, blocks = "blk", statistic = "quadratic",
                    nresample = 199, seed = 5)
  expect_equal(a$p.value, b$p.value)
})

test_that("a seed makes riposte_test reproducible", {
  d <- make_df()
  a <- riposte_test(Y ~ trt | blk, d, statistic = "cauchy", nresample = 199, seed = 99)
  b <- riposte_test(Y ~ trt | blk, d, statistic = "cauchy", nresample = 199, seed = 99)
  expect_equal(a$p.value, b$p.value)
})

test_that("the screen default shrinks and reports lambda and the condition number", {
  d <- make_df()
  r <- riposte_test(Y ~ trt | blk, d, statistic = "screen", nresample = 199, seed = 3)
  expect_equal(r$requested, "screen")
  expect_true(is.finite(r$condition))
  expect_false(is.null(r$screen))
  expect_equal(r$screen$method, "shrink")
  expect_true(r$screen$lambda >= 0 && r$screen$lambda <= 1)
})

test_that("the threshold screen dispatches to the right combination through riposte_test", {
  d <- make_df(effect = 1)
  ## a very low threshold forces the Cauchy; a very high one forces the quadratic
  lo <- riposte_test(Y ~ trt | blk, d, statistic = "screen",
                     screen = riposte_screen_control(method = "threshold", threshold = 1),
                     nresample = 199, seed = 3)
  hi <- riposte_test(Y ~ trt | blk, d, statistic = "screen",
                     screen = riposte_screen_control(method = "threshold", threshold = 1e6),
                     nresample = 199, seed = 3)
  expect_equal(lo$combination, "cauchy")
  expect_equal(hi$combination, "quadratic")
  ## the dispatched p-value matches a direct call to the chosen combination on the
  ## same draws (same RNG kind + seed -> same re-randomizations)
  RNGkind("L'Ecuyer-CMRG"); set.seed(3)
  des <- riposte:::riposte_parse_design(Y ~ trt | blk, d)
  G <- riposte_block_draws(des$z, des$block, 199)
  sm <- riposte_score_matrix(des$y, des$block)
  direct <- riposte_cauchy(sm$scores, des$z, des$block, draws = G)$p.value
  expect_equal(lo$p.value, direct)
})

test_that("riposte_components reports the six representations and the condition number", {
  d <- make_df()
  cmp <- riposte_components(Y ~ trt | blk, d, nresample = 199, seed = 7)
  expect_s3_class(cmp, "riposte_components")
  expect_setequal(names(cmp$component_p), names(riposte_reps_default()))
  expect_true(all(cmp$component_p > 0 & cmp$component_p <= 1))
  expect_true(is.finite(cmp$condition))
  expect_output(print(cmp), "condition number")
})

test_that("a user-supplied representation set flows through riposte_test", {
  d <- make_df(effect = 2)
  reps <- list(raw = function(y) y, abs_dev = function(y) abs(y - median(y)))
  r <- riposte_test(Y ~ trt | blk, d, statistic = "cauchy",
                    representations = reps, nresample = 199, seed = 2)
  expect_s3_class(r, "riposte_test")
  expect_setequal(c(r$kept, r$dropped), c("raw", "abs_dev"))
})
