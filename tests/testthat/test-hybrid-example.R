## The Hybrid is a use of the Cauchy combination: each of the six individual
## tests and the quadratic test contributes one seventh of the statistic.
## Combining the already-combined Cauchy p-value with the quadratic p-value
## at equal weights would instead give the quadratic half the weight.
test_that("the documented Hybrid recipe matches seven coin-derived inputs", {
  skip_if_not_installed("coin")
  d <- asymptotic_example()
  ref <- asymptotic_coin_reference(d)
  individual <- pchisq(as.numeric(coin::statistic(ref, "standardized"))^2,
                       df = 1, lower.tail = FALSE)
  quadratic <- as.numeric(coin::pvalue(ref))
  expected <- pcauchy(sum(tan(pi * (0.5 - c(individual, quadratic)))) / 7,
                      lower.tail = FALSE)

  six <- riposte_test(Y ~ trt | blk, d, statistic = "cauchy",
                      engine = "asymptotic")
  quad <- riposte_test(Y ~ trt | blk, d, statistic = "quadratic",
                       engine = "asymptotic")
  hybrid_inputs <- c(six$component_p, quadratic = quad$p.value)
  hybrid_statistic <- mean(riposte_acat_term(hybrid_inputs))
  hybrid_p <- pcauchy(hybrid_statistic, lower.tail = FALSE)

  expect_length(hybrid_inputs, 7L)
  expect_equal(hybrid_p, expected, tolerance = 1e-9)
  wrong_weights <- pcauchy(mean(riposte_acat_term(c(six$p.value, quad$p.value))),
                           lower.tail = FALSE)
  expect_gt(abs(hybrid_p - wrong_weights), 0.01)
})
