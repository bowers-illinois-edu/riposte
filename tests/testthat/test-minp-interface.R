## test-minp-interface.R
##
## THE POINT. Jake decided on 2026-11-16 (MEMO_paper_changes_2026-11-16.md)
## that the min-p becomes riposte_test()'s default statistic, as the paper
## recommends, and that riposte_minp_control() sets two things about it:
##   pvalue     how each p-value inside the min-p is computed: "midp" (the
##              default until the paper's simulation study of the three rules
##              decides), "randomization", or "asymptotic";
##   quadratic  whether the quadratic combination is one of the p-values. NULL,
##              the default, means yes for a two-sided test and no for a
##              one-sided test, because the quadratic combination responds to
##              departures in both directions. Asking for it with a one-sided
##              alternative is an error whose message says why and names both
##              ways out.
## Under engine = "saddlepoint", which has no min-p, the screen stays the
## default. These tests check the defaults, the errors, covariance adjustment,
## and the printed result. test-minp-rules.R checks the arithmetic.

## ---- riposte_minp_control() ----------------------------------------------------

test_that("riposte_minp_control() defaults to mid-p values and leaves the quadratic to the alternative", {
  ctl <- riposte_minp_control()
  expect_s3_class(ctl, "riposte_minp_control")
  expect_equal(ctl$pvalue, "midp")
  expect_null(ctl$quadratic)
  expect_equal(riposte_minp_control(pvalue = "randomization")$pvalue, "randomization")
  expect_equal(riposte_minp_control(pvalue = "asymptotic")$pvalue, "asymptotic")
  expect_true(riposte_minp_control(quadratic = TRUE)$quadratic)
})

test_that("riposte_minp_control() and riposte_test() reject settings they cannot use", {
  expect_error(riposte_minp_control(pvalue = "exact"))
  expect_error(riposte_minp_control(quadratic = "yes"), "quadratic")
  expect_error(riposte_minp_control(quadratic = NA), "quadratic")
  dat <- minp_interface_design(1)
  expect_error(riposte_test(y ~ z | b, dat, minp = list(pvalue = "midp"), nresample = 49),
               "riposte_minp_control")
})

## ---- the defaults ----------------------------------------------------------------

test_that("the min-p over seven p-values is riposte_test()'s default", {
  dat <- minp_interface_design(2)
  res <- riposte_test(y ~ z | b, dat, nresample = 199, seed = 2)
  expect_equal(res$requested, "minp")
  expect_equal(res$combination, "minp")
  expect_equal(length(res$component_p), 7L)
  expect_equal(names(res$component_p)[7], "quadratic")
  expect_equal(res$minp$pvalue, "midp")
  expect_true(res$minp$quadratic)
  ## the statistic is the smallest of the seven, and the p-value is never
  ## below it (test-minp.R gives the reason)
  expect_equal(res$statistic, min(res$component_p))
  expect_gte(res$p.value, res$statistic - 1e-12)
})

test_that("quadratic = FALSE gives back the six-p-value min-p on the same draws", {
  dat <- minp_interface_design(3)
  res <- riposte_test(y ~ z | b, dat, statistic = "minp", nresample = 199, seed = 3,
                      minp = riposte_minp_control(quadratic = FALSE))
  expect_equal(length(res$component_p), 6L)
  expect_false(res$minp$quadratic)
  ## the same draws riposte_test() makes with this seed
  RNGkind("L'Ecuyer-CMRG"); set.seed(3)
  block <- factor(dat$b)
  G <- riposte_block_draws(dat$z, block, 199L)
  sm <- riposte_score_matrix(dat$y, block)
  old <- riposte:::riposte_minp(sm$scores, dat$z, block, draws = G)
  expect_equal(res$p.value, old$p.value)
})

test_that("a one-sided call uses the min-p without the quadratic combination", {
  dat <- minp_interface_design(4, effect = 1.5)
  g <- riposte_test(y ~ z | b, dat, alternative = "greater", nresample = 199, seed = 4)
  expect_equal(g$requested, "minp")
  expect_equal(length(g$component_p), 6L)
  expect_false("quadratic" %in% names(g$component_p))
  expect_false(g$minp$quadratic)
  ## the effect raises treated outcomes, so "greater" finds it and "less" does not
  l <- riposte_test(y ~ z | b, dat, alternative = "less", nresample = 199, seed = 4)
  expect_lt(g$p.value, l$p.value)
})

test_that("asking for the quadratic combination in a one-sided min-p stops with a message that says why", {
  dat <- minp_interface_design(5)
  msg <- tryCatch(
    riposte_test(y ~ z | b, dat, alternative = "greater", nresample = 49,
                 minp = riposte_minp_control(quadratic = TRUE)),
    error = conditionMessage)
  expect_match(msg, "both directions", fixed = TRUE)
  expect_match(msg, "treated sums fell", fixed = TRUE)
  expect_match(msg, "quadratic = FALSE", fixed = TRUE)
  expect_match(msg, "alternative = \"two.sided\"", fixed = TRUE)
  msg_less <- tryCatch(
    riposte_test(y ~ z | b, dat, alternative = "less", nresample = 49,
                 minp = riposte_minp_control(quadratic = TRUE)),
    error = conditionMessage)
  expect_match(msg_less, "treated sums rose", fixed = TRUE)
})

test_that("each rule runs through riposte_test(), and the randomization rule keeps its guarantee", {
  dat <- minp_interface_design(6)
  block <- factor(dat$b)
  for (rule in c("midp", "randomization", "asymptotic")) {
    res <- riposte_test(y ~ z | b, dat, nresample = 199, seed = 6,
                        minp = riposte_minp_control(pvalue = rule))
    expect_equal(res$minp$pvalue, rule)
    expect_true(res$p.value > 0 && res$p.value <= 1)
    if (rule == "randomization")
      expect_lte(res$p.value, 7 * min(res$component_p) + 1e-12)
    if (rule == "asymptotic") {
      ## the seven p-values are the normal and chi-square tails at the
      ## observed assignment, from the closed-form moments
      sm <- riposte_score_matrix(dat$y, block)
      mom <- riposte_sw_moments(sm$scores, dat$z, block)
      t_obs <- colSums(sm$scores * dat$z)
      zstat <- (t_obs - mom$mu) / sqrt(diag(mom$Sigma))
      cen <- t_obs - mom$mu
      q <- drop(cen %*% solve(mom$Sigma) %*% cen)
      expected <- c(stats::pchisq(zstat^2, 1, lower.tail = FALSE),
                    quadratic = stats::pchisq(q, ncol(sm$scores), lower.tail = FALSE))
      expect_equal(res$component_p, expected, tolerance = 1e-8)
    }
  }
})

## ---- covariance adjustment -------------------------------------------------------

test_that("the min-p's settings carry through covariance adjustment", {
  dat <- minp_interface_design(7)
  two <- riposte_test(y ~ z | b, dat, adjust = ~ x, nresample = 99, seed = 7)
  expect_equal(two$combination, "minp")
  expect_true(two$adjusted)
  expect_equal(length(two$component_p), 7L)
  one <- riposte_test(y ~ z | b, dat, adjust = ~ x, alternative = "greater",
                      nresample = 99, seed = 7)
  expect_equal(length(one$component_p), 6L)
  rand <- riposte_test(y ~ z | b, dat, adjust = ~ x, nresample = 99, seed = 7,
                       minp = riposte_minp_control(pvalue = "randomization"))
  expect_lte(rand$p.value, 7 * min(rand$component_p) + 1e-12)
})

## ---- the other engines -----------------------------------------------------------

test_that("the saddlepoint engine keeps the screen as its default", {
  skip_if_not_installed("fastperm")
  ## two representations, because the quadratic saddlepoint works only for a
  ## few (see test-saddlepoint-engine.R)
  dat <- minp_interface_design(8)
  res <- riposte_test(y ~ z | b, dat, engine = "saddlepoint",
                      representations = list(raw = function(v) v, rank = riposte_rank))
  expect_equal(res$requested, "screen")
})

test_that("the large-sample engine runs the min-p by default and refuses only the screen", {
  ## since 0.0.0.9016 the large-sample min-p draws its reference from the joint
  ## normal approximation (test-minp-normal.R), so the default needs no name
  dat <- minp_interface_design(9)
  res <- riposte_test(y ~ z | b, dat, engine = "asymptotic", nresample = 99, seed = 9)
  expect_equal(res$requested, "minp")
  expect_equal(res$minp$pvalue, "asymptotic")
  msg <- tryCatch(riposte_test(y ~ z | b, dat, statistic = "screen", engine = "asymptotic"),
                  error = conditionMessage)
  expect_match(msg, "screen", fixed = TRUE)
  expect_match(msg, "\"quadratic\"", fixed = TRUE)
  expect_false(grepl("is the max", msg, fixed = TRUE))
})

## ---- the printed result ----------------------------------------------------------

test_that("the printed result names the rule and whether the quadratic combination is included", {
  dat <- minp_interface_design(10)
  res <- riposte_test(y ~ z | b, dat, nresample = 99, seed = 10)
  expect_output(print(res), "mid-p")
  expect_output(print(res), "quadratic combination included")
  one <- riposte_test(y ~ z | b, dat, alternative = "greater", nresample = 99, seed = 10,
                      minp = riposte_minp_control(pvalue = "randomization"))
  expect_output(print(one), "randomization")
  expect_output(print(one), "quadratic combination left out")
})
