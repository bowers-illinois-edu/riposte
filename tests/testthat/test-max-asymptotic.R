## test-max-asymptotic.R
##
## THE POINT. In a large experiment the standardized statistics of the
## representations are jointly normal under the sharp null hypothesis, with the
## null correlation matrix that the design and the outcomes fix. The chance that
## the largest of their absolute values reaches the observed largest can then be
## computed from that multivariate normal distribution, without re-randomizing.
## Rosenbaum (2012, Biometrika) refers his max-T statistic to it, and the paper's
## footnote in Section 4.4 says the software offers it. Jake decided on
## 2026-11-16 that riposte_test(statistic = "max", engine = "asymptotic") gives
## it in every design, not only for matched pairs. These tests check it against
## the one-representation normal test, against a count of normal draws, and
## against the permutation max in a large experiment.

max_normal_design <- function(seed, B = 10L, nb = 40L) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  block <- factor(rep(seq_len(B), each = nb))
  z <- integer(B * nb)
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L
  }
  data.frame(y = rnorm(B * nb, sd = ifelse(z == 1, 1.15, 1)), z = z, b = block)
}

test_that("with one representation the large-sample max is the two-sided normal test", {
  d <- max_normal_design(1)
  raw <- list(raw = function(v) v)
  res <- riposte_test(y ~ z | b, d, statistic = "max", engine = "asymptotic",
                      representations = raw)
  block <- factor(d$b)
  sm <- riposte_score_matrix(d$y, block, raw)
  mom <- riposte_sw_moments(sm$scores, d$z, block)
  zstat <- unname((sum(sm$scores * d$z) - mom$mu) / sqrt(drop(mom$Sigma)))
  expect_equal(res$statistic, abs(zstat), tolerance = 1e-10)
  expect_equal(res$p.value, 2 * stats::pnorm(-abs(zstat)), tolerance = 1e-6)
})

test_that("the large-sample max is the normal chance that the largest |Z| reaches the observed one", {
  skip_if_not_installed("mvtnorm")
  d <- max_normal_design(2)
  res <- riposte_test(y ~ z | b, d, statistic = "max", engine = "asymptotic")
  block <- factor(d$b)
  sm <- riposte_score_matrix(d$y, block)
  mom <- riposte_sw_moments(sm$scores, d$z, block)
  R <- stats::cov2cor(mom$Sigma)
  zstat <- (colSums(sm$scores * d$z) - mom$mu) / sqrt(diag(mom$Sigma))
  expect_equal(res$statistic, max(abs(zstat)), tolerance = 1e-10)
  ## 200,000 draws of the six standardized statistics from N(0, R)
  set.seed(3)
  Zdraw <- matrix(stats::rnorm(2e5 * ncol(R)), ncol = ncol(R)) %*% chol(R)
  share <- mean(apply(abs(Zdraw), 1, max) >= max(abs(zstat)))
  expect_lt(abs(res$p.value - share), 4 * sqrt(share * (1 - share) / 2e5) + 1e-4)
  expect_output(print(res), "approximate")
})

test_that("one-sided large-sample max uses the largest statistic in the named direction", {
  skip_if_not_installed("mvtnorm")
  d <- max_normal_design(4)
  block <- factor(d$b)
  sm <- riposte_score_matrix(d$y, block)
  mom <- riposte_sw_moments(sm$scores, d$z, block)
  R <- stats::cov2cor(mom$Sigma)
  zstat <- (colSums(sm$scores * d$z) - mom$mu) / sqrt(diag(mom$Sigma))
  g <- riposte_test(y ~ z | b, d, statistic = "max", engine = "asymptotic",
                    alternative = "greater")
  expect_equal(g$statistic, max(zstat), tolerance = 1e-10)
  set.seed(5)
  Zdraw <- matrix(stats::rnorm(2e5 * ncol(R)), ncol = ncol(R)) %*% chol(R)
  share <- mean(apply(Zdraw, 1, max) >= max(zstat))
  expect_lt(abs(g$p.value - share), 4 * sqrt(share * (1 - share) / 2e5) + 1e-4)
})

test_that("in a large experiment the large-sample max is close to the permutation max", {
  skip_if_not_installed("mvtnorm")
  d <- max_normal_design(6, B = 20L, nb = 50L)
  asy <- riposte_test(y ~ z | b, d, statistic = "max", engine = "asymptotic")
  per <- riposte_test(y ~ z | b, d, statistic = "max", nresample = 4999, seed = 6)
  expect_lt(abs(asy$p.value - per$p.value), 0.03)
})

test_that("statistics that copy each other count once", {
  skip_if_not_installed("mvtnorm")
  ## the raw outcome twice: the two standardized statistics are equal, their
  ## correlation is 1, and the largest |Z| is that of one representation
  d <- max_normal_design(7)
  twice <- list(a = function(v) v, b = function(v) 2 * v)
  res <- riposte_test(y ~ z | b, d, statistic = "max", engine = "asymptotic",
                      representations = twice)
  once <- riposte_test(y ~ z | b, d, statistic = "max", engine = "asymptotic",
                       representations = list(a = function(v) v))
  expect_equal(res$p.value, once$p.value, tolerance = 1e-6)
})
