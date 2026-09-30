## test-rerandomize.R
##
## The statistical content: the re-randomization must imitate the design exactly
## --- each block keeps its observed treated count in every draw --- and the
## observed assignment must be column 1 so it is exchangeable with the draws.
## Reproducibility comes from the caller's seed.

test_that("column 1 is the observed assignment", {
  set.seed(1)
  z <- rep(c(1, 0), times = 10)
  block <- factor(rep(1:4, each = 5))
  G <- riposte_block_draws(z, block, nresample = 50)
  expect_equal(G[, 1], z, ignore_attr = TRUE)
  expect_equal(ncol(G), 51L)
  expect_equal(nrow(G), length(z))
})

test_that("every draw preserves each block's treated count", {
  set.seed(2)
  z <- c(rep(1, 3), rep(0, 7), rep(1, 2), rep(0, 8))   # blocks with 3/10 and 2/10
  block <- factor(rep(1:2, each = 10))
  G <- riposte_block_draws(z, block, nresample = 200)
  obs_counts <- tapply(z, block, sum)
  ## treated count per block must match the observed count in EVERY column
  for (j in seq_len(ncol(G))) {
    expect_equal(tapply(G[, j], block, sum), obs_counts)
  }
  ## entries are 0/1 only
  expect_true(all(G %in% c(0, 1)))
})

test_that("draws are reproducible under a fixed seed and vary without one", {
  z <- rep(c(1, 0), times = 8)
  block <- factor(rep(1:4, each = 4))
  set.seed(42); G1 <- riposte_block_draws(z, block, nresample = 30)
  set.seed(42); G2 <- riposte_block_draws(z, block, nresample = 30)
  expect_identical(G1, G2)
  ## the draws are not all identical to the observed column (with overwhelming prob)
  expect_false(all(G1[, -1] == z))
})

test_that("nresample = 0 returns just the observed column", {
  z <- c(1, 0, 1, 0)
  block <- factor(c(1, 1, 2, 2))
  G <- riposte_block_draws(z, block, nresample = 0)
  expect_equal(ncol(G), 1L)
  expect_equal(G[, 1], z, ignore_attr = TRUE)
})
