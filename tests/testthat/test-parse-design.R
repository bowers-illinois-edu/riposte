## test-parse-design.R
##
## The two spellings of the design must reduce to the same y, z, block, and the
## treatment must end up coded 0/1 regardless of how the user spelled it.
## riposte_parse_design is internal, so reach it via ::: .

pd <- riposte:::riposte_parse_design

test_that("the formula `| block` spelling and the `blocks` arg agree", {
  set.seed(1)
  d <- data.frame(Y = rnorm(12), trt = rep(c(1, 0), 6),
                  blk = factor(rep(1:3, each = 4)))
  a <- pd(Y ~ trt | blk, d)
  b <- pd(Y ~ trt, d, blocks = "blk")
  expect_equal(a$y, b$y)
  expect_equal(a$z, b$z)
  expect_equal(a$block, b$block)
})

test_that("a two-level factor treatment is coded 0/1 on its second level", {
  d <- data.frame(Y = rnorm(8),
                  arm = factor(rep(c("control", "treated"), 4)),
                  blk = factor(rep(1:2, each = 4)))
  out <- pd(Y ~ arm | blk, d)
  expect_setequal(unique(out$z), c(0, 1))
  ## "treated" sorts after "control", so it becomes the 1 level
  expect_equal(out$z, as.numeric(d$arm == "treated"))
})

test_that("giving the block both ways is an error", {
  d <- data.frame(Y = rnorm(8), trt = rep(c(1, 0), 4),
                  blk = factor(rep(1:2, each = 4)))
  expect_error(pd(Y ~ trt | blk, d, blocks = "blk"), "not both")
})

test_that("a non-binary numeric treatment is rejected", {
  d <- data.frame(Y = rnorm(6), trt = c(0, 1, 2, 0, 1, 2),
                  blk = factor(rep(1:2, each = 3)))
  expect_error(pd(Y ~ trt | blk, d), "0/1")
})
