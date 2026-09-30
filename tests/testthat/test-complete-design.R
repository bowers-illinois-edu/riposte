## Complete randomization fixes the total treated count. Blocking fixes the
## treated count separately in each block. Omitting a block from the public
## formula should mean the first design, equivalent to one explicit block.

test_that("an unblocked formula represents one completely randomized group", {
  d <- asymptotic_example()
  des <- riposte:::riposte_parse_design(Y ~ trt, d)
  expect_equal(des$y, d$Y)
  expect_equal(des$z, d$trt)
  expect_equal(nlevels(des$block), 1L)
  expect_length(des$block, nrow(d))
})

test_that("complete randomization fixes the total but not subgroup counts", {
  d <- asymptotic_example()
  des <- riposte:::riposte_parse_design(Y ~ trt, d)
  set.seed(7)
  draws <- riposte_block_draws(des$z, des$block, nresample = 99)
  expect_equal(draws[, 1], d$trt)
  expect_true(all(colSums(draws) == sum(d$trt)))
  ## d$blk is not in the formula. Its groups must not constrain assignment.
  first_group <- d$blk == levels(d$blk)[1]
  expect_gt(length(unique(colSums(draws[first_group, , drop = FALSE]))), 1L)
})

test_that("blocked randomization fixes each block's treated count", {
  d <- asymptotic_example()
  des <- riposte:::riposte_parse_design(Y ~ trt | blk, d)
  set.seed(7)
  draws <- riposte_block_draws(des$z, des$block, nresample = 99)
  for (ix in split(seq_len(nrow(d)), d$blk)) {
    expect_true(all(colSums(draws[ix, , drop = FALSE]) == sum(d$trt[ix])))
  }
})

test_that("both engines treat no block and one explicit block identically", {
  d <- asymptotic_example()
  d$one <- factor(rep(1L, nrow(d)))
  for (engine in c("permute", "asymptotic")) {
    for (st in c("quadratic", "cauchy")) {
      implicit <- riposte_test(Y ~ trt, d, statistic = st, engine = engine,
                               nresample = 99, seed = 17)
      explicit <- riposte_test(Y ~ trt | one, d, statistic = st, engine = engine,
                               nresample = 99, seed = 17)
      expect_equal(implicit$nblocks, 1L)
      expect_equal(implicit$statistic, explicit$statistic)
      expect_equal(implicit$p.value, explicit$p.value)
    }
  }
})

test_that("both block spellings agree for large-sample calculations", {
  d <- asymptotic_example()
  for (st in c("quadratic", "cauchy")) {
    formula <- riposte_test(Y ~ trt | blk, d, statistic = st,
                            engine = "asymptotic")
    argument <- riposte_test(Y ~ trt, d, blocks = "blk", statistic = st,
                             engine = "asymptotic")
    expect_equal(formula$statistic, argument$statistic)
    expect_equal(formula$p.value, argument$p.value)
  }
})

test_that("component diagnostics also accept complete randomization", {
  d <- asymptotic_example()
  d$one <- factor(rep(1L, nrow(d)))
  implicit <- riposte_components(Y ~ trt, d, nresample = 99, seed = 5)
  explicit <- riposte_components(Y ~ trt | one, d, nresample = 99, seed = 5)
  expect_equal(implicit$nblocks, 1L)
  expect_equal(implicit$Sigma, explicit$Sigma)
  expect_equal(implicit$component_p, explicit$component_p)
})
