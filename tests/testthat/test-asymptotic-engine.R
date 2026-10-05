## The approximation changes the reference used to calculate a p-value, not
## the assignment design or the unadjusted quadratic statistic. Outcomes and
## each block's treated count are held fixed throughout these comparisons.
##
## These tests check arithmetic against coin and hand calculations. They do
## not claim that the approximations have exact size in a finite experiment.

test_that("the quadratic approximation matches coin for unequal blocks", {
  skip_if_not_installed("coin")
  d <- asymptotic_example()
  ref <- asymptotic_coin_reference(d)
  res <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                      engine = "asymptotic")
  expect_equal(res$statistic, as.numeric(coin::statistic(ref, "test")),
               tolerance = 1e-9)
  expect_equal(res$p.value, as.numeric(coin::pvalue(ref)), tolerance = 1e-9)
  expect_equal(res$df, 6L)
})

test_that("the quadratic approximation also matches coin in one block", {
  skip_if_not_installed("coin")
  d <- asymptotic_example()
  d$blk <- factor(rep(1L, nrow(d)))
  ref <- asymptotic_coin_reference(d)
  res <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                      engine = "asymptotic")
  expect_equal(res$statistic, as.numeric(coin::statistic(ref, "test")),
               tolerance = 1e-9)
  expect_equal(res$p.value, as.numeric(coin::pvalue(ref)), tolerance = 1e-9)
})

test_that("one raw score gives the hand-calculated chi-square test", {
  d <- data.frame(Y = c(1, 2, 4, 5, 7, 11),
                  trt = c(0, 0, 0, 1, 1, 1), blk = factor(rep(1, 6)))
  ## The outcome mean is 5. The treated centered sum is 0 + 2 + 6 = 8.
  ## Over assignments of three of six units, its variance is
  ## 3 * 3 / (6 * 5) * (16 + 9 + 1 + 0 + 4 + 36) = 19.8.
  q <- 8^2 / 19.8
  res <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                      representations = list(raw = identity),
                      engine = "asymptotic")
  expect_equal(res$statistic, q)
  expect_equal(res$df, 1L)
  expect_equal(res$p.value, pchisq(q, df = 1, lower.tail = FALSE))
})

test_that("duplicating a score does not add chi-square degrees of freedom", {
  d <- asymptotic_example()
  raw <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                      representations = list(raw = identity),
                      engine = "asymptotic")
  dup <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                      representations = list(raw = identity,
                                             twice = function(y) 2 * y),
                      engine = "asymptotic")
  ## There is still only one independent treated-score sum.
  expect_equal(dup$df, 1L)
  expect_equal(dup$statistic, raw$statistic)
  expect_equal(dup$p.value, raw$p.value)
})

test_that("constant scores do not change the quadratic approximation", {
  d <- asymptotic_example()
  raw <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                      representations = list(raw = identity),
                      engine = "asymptotic")
  res <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                      representations = list(raw = identity,
                                             constant = function(y) rep(1, length(y))),
                      engine = "asymptotic")
  expect_equal(res$dropped, "constant")
  expect_equal(res$df, raw$df)
  expect_equal(res$p.value, raw$p.value)
})

test_that("Cauchy uses the six individual chi-square p-values", {
  skip_if_not_installed("coin")
  d <- asymptotic_example()
  ref <- asymptotic_coin_reference(d)
  z <- as.numeric(coin::statistic(ref, "standardized"))
  p <- pchisq(z^2, df = 1, lower.tail = FALSE)
  ## These p-values are away from zero and one, so the defining tangent
  ## expression can be evaluated directly as an independent reference.
  tc <- mean(tan((0.5 - p) * pi))
  ## cauchy_truncation = 1 is Liu and Xie's original combination, the one this
  ## reference computes; the default truncated version is tested in
  ## test-truncated-cauchy.R
  res <- riposte_test(Y ~ trt | blk, d, statistic = "cauchy",
                      engine = "asymptotic", cauchy_truncation = 1)
  expect_equal(res$statistic, tc, tolerance = 1e-9)
  expect_equal(res$p.value, pcauchy(tc, lower.tail = FALSE), tolerance = 1e-9)
})

test_that("Cauchy on one score reduces to that score's chi-square p-value", {
  d <- asymptotic_example()
  q <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                    representations = list(raw = identity),
                    engine = "asymptotic")
  cauchy <- riposte_test(Y ~ trt | blk, d, statistic = "cauchy",
                         representations = list(raw = identity),
                         engine = "asymptotic")
  expect_equal(cauchy$p.value, q$p.value, tolerance = 1e-12)
})

test_that("small asymptotic p-values are not rounded to zero", {
  ## With 200 zeros assigned control and 200 ones assigned treatment,
  ## the standardized raw statistic squared is N - 1 = 399.
  ## Its chi-square tail is much smaller than machine precision near 1.
  d <- data.frame(Y = rep(0:1, each = 200), trt = rep(0:1, each = 200),
                  blk = factor(rep(1, 400)))
  expected <- pchisq(399, df = 1, lower.tail = FALSE)
  for (st in c("quadratic", "cauchy")) {
    res <- riposte_test(Y ~ trt | blk, d, statistic = st,
                        representations = list(raw = identity),
                        engine = "asymptotic")
    expect_gt(res$p.value, 0)
    ## A relative comparison detects loss of precision in this tiny tail.
    expect_equal(res$p.value / expected, 1, tolerance = 1e-9)
  }
})

test_that("changing the reference leaves the quadratic statistic unchanged", {
  d <- asymptotic_example()
  approx <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                         engine = "asymptotic")
  perm <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                       engine = "permute", nresample = 99, seed = 47)
  expect_equal(approx$statistic, perm$statistic)
  expect_equal(approx$df, perm$df)
})

test_that("ordinary asymptotic calculations do not use random assignments", {
  d <- asymptotic_example()
  set.seed(91)
  before <- .Random.seed
  for (st in c("quadratic", "cauchy")) {
    a <- riposte_test(Y ~ trt | blk, d, statistic = st,
                      engine = "asymptotic", nresample = 19)
    b <- riposte_test(Y ~ trt | blk, d, statistic = st,
                      engine = "asymptotic", nresample = 199)
    expect_identical(a$p.value, b$p.value)
    expect_identical(.Random.seed, before)
  }
})

test_that("permutation remains the default calculation", {
  d <- asymptotic_example()
  for (st in c("quadratic", "cauchy")) {
    default <- riposte_test(Y ~ trt | blk, d, statistic = st,
                            nresample = 99, seed = 17)
    explicit <- riposte_test(Y ~ trt | blk, d, statistic = st,
                             engine = "permute", nresample = 99, seed = 17)
    expect_equal(default$engine, "permute")
    expect_identical(default$p.value, explicit$p.value)
  }
})

test_that("printed asymptotic results identify the approximation", {
  d <- asymptotic_example()
  for (st in c("quadratic", "cauchy")) {
    res <- riposte_test(Y ~ trt | blk, d, statistic = st,
                        engine = "asymptotic")
    expect_equal(res$engine, "asymptotic")
    lines <- capture.output(print(res))
    expect_true(any(grepl("asymptotic|large.sample", lines)))
    expect_false(any(grepl("1999.*re-randomizations", lines)))
  }
})

test_that("cluster approximations use cluster outcomes and assignments", {
  clusters <- asymptotic_example()
  units <- clusters[rep(seq_len(nrow(clusters)), each = 2), ]
  units$cluster <- rep(seq_len(nrow(clusters)), each = 2)
  units$Y <- units$Y + rep(c(-1, 1), nrow(clusters))
  ## The two unit outcomes average to the original cluster outcome. There
  ## must be one randomized observation per cluster, not per person.
  for (st in c("quadratic", "cauchy")) {
    collapsed <- riposte_test(Y ~ trt | blk, clusters, statistic = st,
                              engine = "asymptotic")
    res <- riposte_test(Y ~ trt | blk, units, clusters = "cluster",
                        statistic = st, engine = "asymptotic")
    expect_equal(res$n, nrow(clusters))
    expect_equal(res$statistic, collapsed$statistic)
    expect_equal(res$p.value, collapsed$p.value)
  }
})
