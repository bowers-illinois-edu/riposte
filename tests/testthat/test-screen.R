## test-screen.R
##
## The statistical content:
##   - The threshold screen picks the quadratic when Sigma is well conditioned and
##     the Cauchy when it is not.
##   - The shrink screen leaves a well-conditioned Sigma alone (lambda = 0, the
##     full quadratic), shrinks an ill-conditioned one toward the diagonal, and
##     reaches the diagonal (ignoring the cross-representation dependence) at
##     lambda = 1.
##   - ANCILLARITY (the property the screen's validity rests on): the condition
##     number comes from the closed-form Sigma, which depends on the assignment
##     only through each block's treated COUNT, not through which units are
##     treated. So permuting the assignment within blocks leaves Sigma --- and the
##     condition number --- exactly unchanged. We check this exactly, not
##     approximately.

## a well-conditioned and an ill-conditioned correlation structure
well_cond <- function() diag(4) + 0.05            # nearly identity -> low condition
ill_cond <- function() {
  ## four nearly-collinear representations: condition number large
  R <- matrix(0.98, 4, 4); diag(R) <- 1
  D <- diag(c(1, 2, 0.5, 1.5))
  D %*% R %*% D
}

test_that("threshold screen picks quadratic when well conditioned, Cauchy when not", {
  q <- riposte_screen(well_cond(),
                      riposte_screen_control(method = "threshold", threshold = 100))
  expect_equal(q$choice, "quadratic")
  cc <- riposte_screen(ill_cond(),
                       riposte_screen_control(method = "threshold", threshold = 100))
  expect_equal(cc$choice, "cauchy")
})

test_that("shrink leaves a well-conditioned Sigma alone", {
  Sigma <- well_cond()
  s <- riposte_screen(Sigma, riposte_screen_control(method = "shrink", threshold = 100))
  expect_equal(s$lambda, 0)
  expect_equal(s$Sigma, Sigma, ignore_attr = TRUE)
})

test_that("shrink reduces the condition number of an ill-conditioned Sigma", {
  Sigma <- ill_cond()
  s <- riposte_screen(Sigma, riposte_screen_control(method = "shrink", threshold = 50))
  expect_gt(s$lambda, 0)
  expect_lt(s$condition_shrunk, s$condition)
  ## the achieved condition number is at (or below) the target, up to rounding
  expect_lte(s$condition_shrunk, 50 * 1.0001)
})

test_that("shrink limits: lambda = 0 is the full Sigma, lambda = 1 is diagonal", {
  Sigma <- ill_cond()
  full <- riposte_screen(Sigma, riposte_screen_control(method = "shrink", lambda = 0))
  expect_equal(full$Sigma, Sigma, ignore_attr = TRUE)
  diagonal <- riposte_screen(Sigma, riposte_screen_control(method = "shrink", lambda = 1))
  expect_equal(diagonal$Sigma, diag(diag(Sigma)), ignore_attr = TRUE)
})

test_that("the condition number is ancillary: invariant to the within-block assignment", {
  ## a real design; Sigma comes from the closed-form moments
  set.seed(20260623)
  block <- factor(rep(1:5, each = 10))
  y <- rnorm(length(block))
  sm <- riposte_score_matrix(y, block)

  ## several different assignments with the SAME treated count per block
  cond_numbers <- replicate(20, {
    z <- integer(length(block))
    for (b in levels(block)) {
      ix <- which(block == b); z[ix][sample.int(length(ix), 5L)] <- 1L
    }
    Sigma <- riposte_sw_moments(sm$scores, z, block)$Sigma
    riposte_screen(Sigma)$condition
  })
  ## the closed-form Sigma does not depend on WHICH units are treated, only the
  ## per-block count, so every condition number is identical
  expect_lt(stats::sd(cond_numbers), 1e-9)
})
