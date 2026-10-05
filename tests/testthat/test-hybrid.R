## test-hybrid.R
##
## THE POINT. The Hybrid gives equal weight to seven p-values: each of the six
## representations tested alone and the quadratic combination. Its documented
## recipe (test-hybrid-example.R) takes all seven from the large-sample
## approximation, and that version returns exactly 1 whenever one
## representation's treated sum equals its re-randomization mean exactly: that
## representation's p-value is 1 and its Cauchy term tan((0.5 - 1) * pi) is minus
## infinity. On a coarse outcome with a few dozen units this happens (Wallsten and
## Nteta 2016, 36 respondents, in the Sarkar and Coppock reanalysis).
##
## statistic = "hybrid" refers the Hybrid to the randomization distribution, like
## the quadratic, Cauchy, and max: at every re-randomization it averages the
## Cauchy terms of the six representations' two-sided mid-p values and the
## quadratic's mid-p value, and the p-value is the share of re-randomizations
## whose average is at least the observed one. The mid-p value counts ties with
## the observed statistic as half, so it stays below 1.

test_that("the hybrid is the share of re-randomizations whose seven-term Cauchy average is at least the observed one", {
  ## rebuild the statistic by hand from the same draws riposte_test() uses
  d <- wallsten_36(); b <- factor(rep(1, nrow(d)))
  RNGkind("L'Ecuyer-CMRG"); set.seed(1)
  sm <- riposte_score_matrix(d$Y, b)
  mom <- riposte_sw_moments(sm$scores, d$Z, b)
  Tm <- riposte_perm_stats(sm$scores, d$Z, b, draws = riposte_block_draws(d$Z, b, 1999L))$stats
  inv <- riposte_std_pinv(mom$Sigma)
  cs <- sweep(sweep(Tm, 2, mom$mu), 2, inv$sd, "/")
  Q <- rowSums((cs %*% inv$R_pinv) * cs)
  P <- cbind(apply(Tm, 2, riposte_midp), quadratic = riposte_midp(Q, alternative = "greater"))
  h <- rowMeans(riposte_acat_term(P))
  by_hand <- mean(h >= h[1] - riposte_tie_tol(h[1]))

  res <- riposte_test(Y ~ Z, data = d, statistic = "hybrid", seed = 1)
  expect_equal(res$p.value, by_hand)
  expect_equal(res$combination, "hybrid")
  expect_named(res$component_p, c(colnames(Tm), "quadratic"))
})

test_that("a representation whose treated sum sits exactly at its mean does not force the hybrid to 1", {
  d <- wallsten_36()
  expect_equal(mean(pmax(d$Y, 1 - d$Y)[d$Z == 1]), 31 / 36)
  expect_equal(mean(pmax(d$Y, 1 - d$Y)), 31 / 36)
  ## the documented large-sample recipe gives exactly 1
  six <- riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic")
  quad <- riposte_test(Y ~ Z, d, statistic = "quadratic", engine = "asymptotic")
  large_sample <- pcauchy(mean(riposte_acat_term(c(six$component_p, quad$p.value))), lower.tail = FALSE)
  expect_equal(large_sample, 1)
  ## the raw, rank, and Huber p-values are each near 0.02 and the difference in
  ## means gives 0.028; the re-randomization hybrid gives about 0.05 (0.055 with
  ## seed 1 and these rows in this order), not 1
  expect_lt(riposte_test(Y ~ Z, d, statistic = "hybrid", seed = 1)$p.value, 0.1)
})

test_that("on a binary outcome the hybrid is the permutation difference in means", {
  ## every representation of a 0/1 outcome is a + bY, so all seven terms are the
  ## same mid-p value at every re-randomization and the hybrid orders the
  ## re-randomizations as the difference in means does
  RNGkind("L'Ecuyer-CMRG"); set.seed(8)
  d <- data.frame(Z = rep(0:1, 100), Y = rbinom(200, 1, 0.4))
  raw_only <- list(raw = function(y) y)
  expect_equal(riposte_test(Y ~ Z, d, statistic = "hybrid", seed = 1)$p.value,
               riposte_test(Y ~ Z, d, statistic = "quadratic", representations = raw_only,
                            seed = 1)$p.value)
})

test_that("under the sharp null the hybrid holds the 0.05 level", {
  skip_on_cran()
  ## 300 data sets with no effect on a coarse 5-point outcome, 36 units, half
  ## treated; SE of the rejection rate is sqrt(0.05 * 0.95 / 300) = 0.0126
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261004)
  rej <- vapply(seq_len(300), function(i) {
    d <- data.frame(Y = sample(c(0, 0.25, 0.5, 0.75, 1), 36, replace = TRUE),
                    Z = sample(rep(0:1, 18)))
    riposte_test(Y ~ Z, d, statistic = "hybrid", nresample = 199L)$p.value <= 0.05
  }, logical(1))
  expect_lt(mean(rej), 0.05 + 3 * sqrt(0.05 * 0.95 / 300))
})

test_that("the large-sample hybrid uses the truncated Cauchy combination of its seven p-values", {
  ## With engine = "asymptotic" the seven p-values come from the chi-square
  ## approximations and are combined with Gui, Jiang, and Wang's truncated
  ## conversion (see test-truncated-cauchy.R), so a p-value of 1 cannot set the
  ## hybrid to 1.
  d <- wallsten_36()
  six <- riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic")
  quad <- riposte_test(Y ~ Z, d, statistic = "quadratic", engine = "asymptotic")
  res <- riposte_test(Y ~ Z, d, statistic = "hybrid", engine = "asymptotic")
  expect_equal(res$p.value, riposte_truncated_cauchy(c(six$component_p, quadratic = quad$p.value)))
  expect_lt(res$p.value, 0.1)
})

test_that("the hybrid has no saddlepoint version and no one-sided form", {
  d <- wallsten_36()
  expect_error(riposte_test(Y ~ Z, d, statistic = "hybrid", engine = "saddlepoint"), "permute")
  expect_error(riposte_test(Y ~ Z, d, statistic = "hybrid", alternative = "greater"), "one-sided")
})
