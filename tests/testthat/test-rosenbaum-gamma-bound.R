## The upper bound on the p-value of testing one hypothesis twice when gamma > 1.
##
## Under Rosenbaum's sensitivity model with hidden bias at most gamma, the signs of
## the matched-pair differences are independent, and each is positive with some
## chance pi_i <= theta = gamma / (1 + gamma). Every pair score is nonnegative, so
## turning a sign positive can only raise the largest standardized statistic. The
## chance that this maximum reaches its observed value is therefore largest when
## every pi_i = theta, and that chance is the bound. POINT.md asks riposte to compute
## it from the distribution of the signs at theta (every sign pattern listed, or
## patterns drawn at random) and to use the normal approximation only under
## engine = "asymptotic".

make_pairs <- function(I, seed, shift = 0.4) {
  set.seed(seed)
  d <- data.frame(pair = rep(seq_len(I), each = 2), z = rep(0:1, I))
  d$y <- stats::rnorm(2 * I) + shift * d$z
  d
}
pair_d <- function(d) { w <- split(d, d$pair); vapply(w, function(a) a$y[a$z == 1] - a$y[a$z == 0], numeric(1)) }
pair_scores <- function(y, reps = riposte_pair_reps()) sapply(reps, function(f) attr(f, "pair_scores")(y))

## Independent brute force: loop over every sign pattern, give it its chance when
## each sign is positive with chance theta, and add up the chances of the patterns
## whose largest standardized statistic reaches the observed one.
exact_bound <- function(y, H, gamma) {
  theta <- gamma / (1 + gamma)
  keep <- y != 0
  y <- y[keep]; H <- H[keep, colSums(H[keep, , drop = FALSE]) > 0, drop = FALSE]
  I <- length(y)
  mu <- theta * colSums(H)
  sdv <- sqrt(theta * (1 - theta) * colSums(H^2))
  max_z <- function(b) max((as.vector(b %*% H) - mu) / sdv)
  obs <- max_z(as.numeric(y > 0))
  total <- 0
  for (i in 0:(2^I - 1)) {
    b <- as.numeric(bitwAnd(i, 2^(0:(I - 1))) > 0)
    if (max_z(b) >= obs - 1e-12) total <- total + theta^sum(b) * (1 - theta)^(I - sum(b))
  }
  total
}

test_that("listing sign patterns at theta gives the brute-force bound", {
  set.seed(31)
  y <- stats::rnorm(10) + 0.5
  H <- pair_scores(y)
  for (g in c(1, 1.5, 3))
    expect_equal(riposte_pairs_perm_test(y, H, nresample = NULL, gamma = g)$p.value,
                 exact_bound(y, H, g), tolerance = 1e-12, label = paste("gamma", g))
})

test_that("gamma = 1, the default, leaves the sign-flip test as it was", {
  set.seed(33)
  y <- stats::rnorm(9) + 0.3
  H <- pair_scores(y)
  expect_identical(riposte_pairs_perm_test(y, H, nresample = NULL, gamma = 1),
                   riposte_pairs_perm_test(y, H, nresample = NULL))
})

test_that("the bound holds for every pattern of bias no larger than gamma", {
  ## With |y| fixed at 1, ..., 8 the scores do not depend on the signs, so the bound
  ## can be computed for all 2^8 sign patterns. If each pattern's true chance comes
  ## from pi_i <= theta, the chance that the bound is at most a must be at most a.
  I <- 8; gamma <- 2; theta <- gamma / (1 + gamma)
  S <- as.matrix(expand.grid(rep(list(c(-1, 1)), I)))
  H <- pair_scores(seq_len(I))
  p <- apply(S, 1, function(s) riposte_pairs_perm_test(s * seq_len(I), H, nresample = NULL, gamma = gamma)$p.value)
  biases <- list(largest = rep(theta, I), none = rep(0.5, I), mixed = rep(c(theta, 0.5), I / 2))
  for (nm in names(biases)) {
    pos <- biases[[nm]]                       # each pair's chance of a positive sign
    chance <- apply(S, 1, function(s) prod(ifelse(s > 0, pos, 1 - pos)))
    for (a in c(0.05, 0.10, 0.20, 0.50))
      expect_lte(sum(chance[p <= a + 1e-12]), a + 1e-12, label = paste(nm, "bias, level", a))
  }
})

test_that("signs drawn at theta agree with the listed patterns within simulation error", {
  set.seed(34)
  y <- stats::rnorm(14) + 0.4
  H <- pair_scores(y)
  ex <- riposte_pairs_perm_test(y, H, nresample = NULL, gamma = 1.5)$p.value
  set.seed(35)
  mc <- riposte_pairs_perm_test(y, H, nresample = 4000, gamma = 1.5)$p.value
  expect_lt(abs(ex - mc), 4 * sqrt(ex * (1 - ex) / 4000) + 0.005)
  ## the estimate counts the observed signs as one of the draws
  expect_equal(mc * 4001, round(mc * 4001), tolerance = 1e-9)
})

test_that("riposte_test() at gamma > 1 uses the sign distribution unless engine = 'asymptotic'", {
  skip_if_not_installed("mvtnorm")
  d <- make_pairs(10, 32, shift = 0.6)
  dd <- pair_d(d)
  H <- pair_scores(dd)
  fit <- riposte_test(y ~ z | pair, d, gamma = 1.5, nresample = 2^10, seed = 1)
  expect_true(fit$enumerated)
  expect_equal(fit$p.value, exact_bound(dd, H, 1.5), tolerance = 1e-12)
  ## the statistic is the largest deviate standardized at theta, as in the normal approximation
  bnd <- riposte_pairs_bound(dd, H, gamma = 1.5)
  expect_equal(fit$statistic, max(bnd$deviates), tolerance = 1e-12)
  asy <- riposte_test(y ~ z | pair, d, gamma = 1.5, engine = "asymptotic")
  expect_equal(asy$p.value, bnd$p.value)
})

test_that("with many pairs the default engine draws signs at theta, reproducibly", {
  d <- make_pairs(40, 36, shift = 0.3)
  a <- riposte_test(y ~ z | pair, d, gamma = 2, nresample = 999, seed = 1)
  b <- riposte_test(y ~ z | pair, d, gamma = 2, nresample = 999, seed = 1)
  expect_false(a$enumerated)
  expect_identical(a$p.value, b$p.value)
  expect_equal(a$p.value * 1000, round(a$p.value * 1000), tolerance = 1e-9)
})

test_that("alternative = 'less' at gamma > 1 mirrors 'greater' on the negated outcome", {
  d <- make_pairs(10, 37)
  d2 <- d; d2$y <- -d2$y
  a <- riposte_test(y ~ z | pair, d,  gamma = 2, nresample = 2^10, alternative = "less")
  b <- riposte_test(y ~ z | pair, d2, gamma = 2, nresample = 2^10, alternative = "greater")
  expect_equal(a$p.value, b$p.value)
})

test_that("the printed result says which calculation gave the bound", {
  skip_if_not_installed("mvtnorm")
  d <- make_pairs(10, 38)
  exact <- paste(utils::capture.output(print(riposte_test(y ~ z | pair, d, gamma = 2, nresample = 2^10))), collapse = "\n")
  approx <- paste(utils::capture.output(print(riposte_test(y ~ z | pair, d, gamma = 2, engine = "asymptotic"))), collapse = "\n")
  expect_match(exact, "upper bound")
  expect_no_match(exact, "normal|approximat")
  expect_match(approx, "approximat")
})
