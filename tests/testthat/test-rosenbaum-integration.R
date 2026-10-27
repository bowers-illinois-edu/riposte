## Rosenbaum (2012) matched-pair tests reached through riposte_test(): the pair
## representations must reproduce, through riposte's own machinery (centred scores,
## Strasser-Weber moments, re-randomization), what the standalone riposte_pairs_*()
## functions and Rosenbaum's published numbers say.

make_pairs <- function(I, seed, shift = 0.4) {
  set.seed(seed)
  d <- data.frame(pair = rep(seq_len(I), each = 2), z = rep(0:1, I))
  d$y <- stats::rnorm(2 * I) + shift * d$z
  d
}
pair_d <- function(d) { w <- split(d, d$pair); vapply(w, function(a) a$y[a$z == 1] - a$y[a$z == 0], numeric(1)) }

test_that("pair representations are Strasser-Weber linear statistics equal to sum B_i h_i - sum h_i / 2", {
  d <- make_pairs(12, 11)
  reps <- riposte_pair_reps()
  sm <- riposte_score_matrix(d$y, d$pair, reps)
  Tobs <- as.vector(crossprod(d$z, sm$scores))           # sum over treated units of centred scores
  dd <- pair_d(d)
  for (j in seq_along(reps)) {
    h <- attr(reps[[j]], "pair_scores")(dd)
    expect_equal(Tobs[j], sum((dd > 0) * h) - sum(h) / 2, tolerance = 1e-12, label = names(reps)[j])
  }
})

test_that("closed-form moments equal Rosenbaum's worst-case moments at Gamma = 1", {
  d <- make_pairs(15, 12)
  reps <- riposte_pair_reps()
  sm <- riposte_score_matrix(d$y, d$pair, reps)
  mom <- riposte_sw_moments(sm$scores, d$z, d$pair)
  dd <- pair_d(d)
  H <- sapply(reps, function(f) attr(f, "pair_scores")(dd))
  expect_equal(unname(mom$mu), rep(0, ncol(H)), tolerance = 1e-12)
  expect_equal(unname(mom$Sigma), unname(crossprod(H) / 4), tolerance = 1e-10)   # theta (1 - theta) = 1/4
})

test_that("riposte_test(gamma = 1) equals the independent sign-flip enumeration", {
  I <- 10
  d <- make_pairs(I, 13, shift = 0.8)
  dd <- pair_d(d)
  reps <- riposte_pair_reps()
  H <- sapply(reps, function(f) attr(f, "pair_scores")(dd))
  ## independent brute force over all 2^I sign patterns
  S <- t(as.matrix(expand.grid(rep(list(0:1), I))))
  Z <- sweep(sweep(crossprod(S, H), 2, colSums(H) / 2), 2, sqrt(colSums(H^2) / 4), "/")
  mx <- apply(Z, 1, max)
  obs <- max((colSums((dd > 0) * H) - colSums(H) / 2) / sqrt(colSums(H^2) / 4))
  p_exact <- mean(mx >= obs - 1e-12)
  fit <- riposte_test(y ~ z | pair, d, statistic = "max", gamma = 1, nresample = 2^I, seed = 1)
  expect_equal(fit$statistic, obs, tolerance = 1e-12)
  expect_equal(fit$p.value, p_exact)
  expect_match(paste(utils::capture.output(print(fit)), collapse = "\n"), "exact")
})

test_that("Monte-Carlo sign flips agree with the exact p-value within simulation error", {
  d <- make_pairs(14, 14, shift = 0.5)
  ex  <- riposte_test(y ~ z | pair, d, statistic = "max", gamma = 1, nresample = 2^14, seed = 1)$p.value
  mc  <- riposte_test(y ~ z | pair, d, statistic = "max", gamma = 1, nresample = 4000, seed = 2)$p.value
  expect_lt(abs(ex - mc), 4 * sqrt(ex * (1 - ex) / 4000) + 0.005)
})

test_that("Gamma > 1 through riposte_test equals riposte_pairs_bound and CRAN testtwice", {
  skip_if_not_installed("testtwice"); skip_if_not_installed("mvtnorm")
  e <- new.env(); utils::data("smokerlead", package = "testtwice", envir = e)
  y <- e$smokerlead$lead
  d <- data.frame(pair = rep(seq_along(y), each = 2), z = rep(0:1, length(y)),
                  y = as.vector(rbind(0, y)))                 # control = 0, treated = the difference
  for (g in c(1, 2, 3)) {
    fit <- riposte_test(y ~ z | pair, d, statistic = "max", gamma = g, engine = "asymptotic",
                        representations = riposte_pair_reps(c("u878", "u867")))
    ref <- testtwice::testtwice(y, gamma = g, u878 = TRUE, u867 = TRUE)
    expect_lte(abs(fit$p.value - ref$pval), 0.005, label = paste("Gamma", g))
    expect_equal(fit$statistic, max(ref$dev), tolerance = 0.02, label = paste("max deviate, Gamma", g))
  }
})

test_that("adaptive test: p-value is at least the smallest rejecting level (brute-force scan)", {
  set.seed(15)
  for (rep in 1:6) {
    I <- 12
    y <- stats::rnorm(I) + 0.7
    ct <- riposte_adaptive_counts(y)
    g <- riposte_adaptive_grid(ct$I1, ct$I2, 1)
    alphas <- sort(unique(g$U[g$U > 0 & g$U < 1]))
    rej <- vapply(alphas, function(a) { cr <- riposte_adaptive_pick(g, a); !is.null(cr) && (ct$B1 >= cr$kB || ct$T >= cr$kT) }, logical(1))
    p_scan <- if (any(rej)) min(alphas[rej]) else 1
    ## equal when the regions are nested; larger when they are not (see the validity test)
    expect_gte(riposte_adaptive_pvalue(y, 1), p_scan - 1e-12)
  }
})

test_that("adaptive test p-value is valid: P(p <= alpha) <= alpha under all 2^I sign patterns", {
  I <- 10
  mag <- seq_len(I)                                 # untied |differences|
  S <- as.matrix(expand.grid(rep(list(c(-1, 1)), I)))
  for (gam in c(1, 1.5)) {
    ## under the sensitivity bound each sign is positive with probability kappa
    kap <- gam / (1 + gam)
    pr <- apply(S, 1, function(s) prod(ifelse(s > 0, kap, 1 - kap)))
    pv <- apply(S, 1, function(s) riposte_adaptive_pvalue(s * mag, gam, scan = Inf))
    for (a in c(0.01, 0.05, 0.10, 0.25))
      expect_lte(sum(pr[pv <= a + 1e-12]), a + 1e-12)
  }
})

test_that("the adaptive test through riposte_test() matches the standalone function", {
  d <- make_pairs(40, 16, shift = 0.6)
  dd <- pair_d(d)
  fit <- riposte_test(y ~ z | pair, d, statistic = "adaptive", gamma = 1.2, alpha = 0.05)
  expect_equal(fit$p.value, riposte_adaptive_pvalue(dd, 1.2))
  expect_equal(fit$adaptive$reject, riposte_adaptive_test(dd, 1.2, 0.05)$reject)
  less <- riposte_test(y ~ z | pair, d, statistic = "adaptive", alternative = "less")
  expect_equal(less$p.value, riposte_adaptive_pvalue(-dd, 1))
})

test_that("alternative = 'less' mirrors 'greater' on the negated outcome", {
  d <- make_pairs(12, 17)
  d2 <- d; d2$y <- -d2$y
  a <- riposte_test(y ~ z | pair, d,  statistic = "max", gamma = 2, engine = "asymptotic", alternative = "less")
  b <- riposte_test(y ~ z | pair, d2, statistic = "max", gamma = 2, engine = "asymptotic", alternative = "greater")
  expect_equal(a$p.value, b$p.value)
})

test_that("matched-pair options reject designs and arguments they do not support", {
  d <- make_pairs(8, 18)
  d3 <- data.frame(y = rnorm(9), z = rep(c(0, 0, 1), 3), b = rep(1:3, each = 3))
  expect_error(riposte_test(y ~ z | b, d3, statistic = "max", gamma = 1, representations = riposte_pair_reps()), "matched pairs")
  expect_error(riposte_test(y ~ z | pair, d, gamma = 0.5), ">= 1")
  expect_error(riposte_test(y ~ z | pair, d, statistic = "quadratic", gamma = 1), "max")
  expect_error(riposte_test(y ~ z | pair, d, statistic = "adaptive", clusters = "pair"), "do not support")
  expect_error(riposte_test(y ~ z | pair, d, statistic = "max", gamma = 2, representations = riposte_reps_default()), "riposte_pair_reps")
})

test_that("pair representations also work in the usual interface (two-sided) and in components", {
  d <- make_pairs(25, 19, shift = 0.7)
  fit <- riposte_test(y ~ z | pair, d, statistic = "quadratic", representations = riposte_pair_reps(),
                      nresample = 500, seed = 1)
  expect_true(fit$p.value > 0 && fit$p.value <= 1)
  comp <- riposte_components(y ~ z | pair, d, representations = riposte_pair_reps(), nresample = 200, seed = 1)
  expect_equal(sort(comp$kept), sort(names(riposte_pair_reps())))
})

test_that("identical statistics (few pairs) are left out of the orthant probability, not an error", {
  skip_if_not_installed("mvtnorm")
  ## with six pairs u878 and u867 give identical scores, so the correlation matrix is singular
  d6 <- c(2, 6, -0.5, 8, 0.5, 4.5)
  H <- do.call(cbind, lapply(riposte_pair_reps(), function(f) attr(f, "pair_scores")(d6)))
  expect_equal(H[, "u878"], H[, "u867"])
  full <- riposte_pairs_bound(d6, H, gamma = 2)
  expect_true(full$p.value > 0 && full$p.value < 1)
  ## dropping the copy by hand gives the same probability
  hand <- riposte_pairs_bound(d6, H[, colnames(H) != "u867"], gamma = 2)
  expect_equal(full$p.value, hand$p.value, tolerance = 1e-8)
  ## and the interface runs on such a design
  d <- data.frame(y = c(5,3, 8,2, 4,4.5, 9,1, 6,5.5, 7,2.5),
                  z = rep(c(1, 0), 6), pair = rep(1:6, each = 2))
  fit <- riposte_test(y ~ z | pair, d, statistic = "max", gamma = 2, engine = "asymptotic")
  expect_equal(fit$p.value, full$p.value, tolerance = 1e-8)
})

test_that("matched-pair test with six pairs runs through riposte_test()", {
  set.seed(1)
  d <- data.frame(y = c(rnorm(6) + 1, rep(0, 6)), trt = rep(c(1, 0), each = 6), pair = rep(1:6, 2))
  r <- riposte_test(y ~ trt | pair, d, statistic = "max", gamma = 1)
  expect_false(is.nan(r$p.value)); expect_true(is.finite(r$p.value))
})

test_that("riposte_test() keeps its old positional signature", {
  fm <- names(formals(riposte_test))
  expect_lt(match("alternative", fm), match("gamma", fm))
  expect_lt(match("alternative", fm), match("alpha", fm))
  expect_equal(fm[seq_len(match("alternative", fm))][1:3], c("formula", "data", "blocks"))
})

