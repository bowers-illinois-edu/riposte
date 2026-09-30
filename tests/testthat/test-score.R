## test-score.R
##
## The statistical content: each representation becomes a within-block-centred
## score, so every column has mean zero within every block (the Strasser-Weber
## condition that makes the permutation mean of the linear statistic zero). A
## representation constant within every block carries no information and must be
## dropped. A representation returning the wrong length is a user error and must
## be caught.

test_that("every kept score column has within-block mean zero", {
  set.seed(1)
  y <- rnorm(30)
  block <- factor(rep(1:3, each = 10))
  sm <- riposte_score_matrix(y, block)
  for (col in seq_len(ncol(sm$scores))) {
    bm <- tapply(sm$scores[, col], block, mean)
    expect_true(all(abs(bm) < 1e-10))
  }
})

test_that("a within-block-constant representation is dropped", {
  set.seed(2)
  y <- rnorm(20)
  block <- factor(rep(1:2, each = 10))
  reps <- c(riposte_reps_default(), list(const = function(v) rep(1, length(v))))
  sm <- riposte_score_matrix(y, block, representations = reps)
  expect_true("const" %in% sm$dropped)
  expect_false("const" %in% sm$kept)
  expect_false("const" %in% colnames(sm$scores))
})

test_that("a representation returning the wrong length is an error", {
  y <- rnorm(10)
  block <- factor(rep(1:2, each = 5))
  bad <- list(bad = function(v) v[-1])   # returns length n-1
  expect_error(riposte_score_matrix(y, block, representations = bad),
               "one score per unit")
})

test_that("scores have one row per unit and one column per kept representation", {
  set.seed(3)
  y <- rnorm(24)
  block <- factor(rep(1:4, each = 6))
  sm <- riposte_score_matrix(y, block)
  expect_equal(nrow(sm$scores), length(y))
  expect_equal(ncol(sm$scores), length(sm$kept))
  expect_setequal(c(sm$kept, sm$dropped), names(riposte_reps_default()))
})
