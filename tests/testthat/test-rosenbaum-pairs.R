## Tests for R/rosenbaum-pairs.R. Targets are numbers reported in Rosenbaum (2012a,
## Biometrika; 2012b, Ann. Appl. Stat.), and an independent brute-force enumeration.

test_that("adaptive test: Annals Section 4.2 worked example (I = 250, Gamma = 4)", {
  y <- seq_len(250) * c(1, -1)               # untied |y| = 1..250
  cr <- riposte_adaptive_critical(I1 = 84, I2 = 83, gamma = 4, alpha = 0.05)
  expect_equal(c(cr$kB, cr$kT), c(74, 216))
  expect_equal(round(cr$size, 4), 0.0488)
  expect_equal(round(cr$pB, 4), 0.0370)
  expect_equal(round(cr$pT, 4), 0.0320)
  out <- riposte_adaptive_test(y, gamma = 4)
  expect_equal(c(out$I1, out$I2), c(84, 83))
})

test_that("adaptive test: Annals Section 5 narrative (Gamma = 1.96, I1 = 110, I2 = 106)", {
  ## The paper prints I1 + I2 = 226; 110 + 106 = 216. The probabilities below are
  ## reproduced with (110, 106), so '226' appears to be a typo for 216.
  cr <- riposte_adaptive_critical(110, 106, gamma = 1.96, alpha = 0.05)
  expect_equal(c(cr$kB, cr$kT), c(82, 237))
  expect_equal(round(cr$pB, 4), 0.0381)
  expect_equal(round(cr$pT, 4), 0.0293)
  expect_equal(round(cr$size, 4), 0.0475)
})

test_that("adaptive test has exact size under full sign enumeration (gamma = 1)", {
  I <- 12
  y <- seq_len(I)                              # |y| = 1..I, all positive for the observed
  ## bands from ranks (I = 12, lambda = 1/3, 2/3): ranks 8..12 and 4..7
  q <- seq_len(I)
  b1 <- q >= (1 - 1/3) * I; b2 <- q < (1 - 1/3) * I & q >= (1 - 2/3) * I
  cr <- riposte_adaptive_critical(I1 = sum(b1), I2 = sum(b2), gamma = 1, alpha = 0.10)
  ## brute force over all 2^12 sign vectors
  S <- as.matrix(expand.grid(rep(list(0:1), I)))
  B1 <- S[, b1] %*% rep(1, sum(b1)); B2 <- S[, b2] %*% rep(1, sum(b2))
  rej <- (B1 >= cr$kB) | (2 * B1 + B2 >= cr$kT)
  expect_equal(mean(rej), cr$size, tolerance = 1e-12)
  expect_lte(cr$size, 0.10)
})

test_that("equicorrelated-normal cutoffs: Biometrika Table 4", {
  skip_if_not_installed("mvtnorm")
  eq <- function(L, rho) { R <- matrix(rho, L, L); diag(R) <- 1; R }
  k <- mvtnorm::qmvnorm(0.95, corr = eq(10, 0.9), tail = "lower.tail")$quantile
  expect_equal(round(k, 3), 2.077, tolerance = 0.003)
  expect_equal(round(10 * (1 - stats::pnorm(k)), 3), 0.189, tolerance = 0.01)
  expect_equal(round(stats::qnorm(1 - 0.05 / 10), 3), 2.576)
})

test_that("scores: shapes, zero handling, scale", {
  y <- c(-3, 0, 1, 5, -2, 4)
  for (ty in c("wilcoxon", "brown", "noether")) {
    h <- riposte_pairs_scores(y, ty)
    expect_length(h, length(y)); expect_equal(h[2], 0); expect_equal(max(h), 1)
  }
  h <- riposte_pairs_scores(y, "ustat", m = 8, m1 = 7, m2 = 8)
  expect_equal(h[2], 0)
  expect_error(riposte_pairs_scores(y, "ustat", m = 8, m1 = 9, m2 = 8))
})

test_that("exact sign-flip test equals independent enumeration at gamma = 1", {
  set.seed(1)
  y <- rnorm(10) + 0.5
  H <- cbind(riposte_pairs_scores(y, "wilcoxon"), riposte_pairs_scores(y, "brown"))
  out <- riposte_pairs_perm_test(y, H, nresample = NULL)
  ## independent: loop over all 2^10 sign vectors
  obs <- NULL; sims <- numeric(2^10)
  mx <- function(b) { T <- b %*% H; z <- (T - colSums(H) / 2) / sqrt(colSums(H^2) / 4); max(z) }
  obs <- mx(as.numeric(y > 0))
  for (i in 0:(2^10 - 1)) sims[i + 1] <- mx(as.numeric(bitwAnd(i, 2^(0:9)) > 0))
  expect_equal(out$p.value, mean(sims >= obs - 1e-12))
})

test_that("matches the CRAN package testtwice on its smokerlead data (Biometrika Tables 2-3)", {
  skip_if_not_installed("testtwice"); skip_if_not_installed("mvtnorm")
  data("smokerlead", package = "testtwice")
  y <- smokerlead$lead
  ## Table 2, Wilcoxon-free: (8,7,8) alone and (8,6,7) alone, normal approximation
  one <- function(h, g) { th <- g / (1 + g); 1 - stats::pnorm((sum((y > 0) * h) - th * sum(h)) / sqrt(th * (1 - th) * sum(h^2))) }
  u878 <- riposte_pairs_scores(y, "ustat", 8, 7, 8); u867 <- riposte_pairs_scores(y, "ustat", 8, 6, 7)
  expect_equal(round(one(u878, 3), 3), 0.236)
  expect_equal(round(one(u867, 3), 3), 0.050)
  ## joint procedure vs testtwice::tt (same normal-orthant calculation)
  ref <- testtwice::testtwice(y, gamma = 3, u878 = TRUE, u867 = TRUE)
  ours <- riposte_pairs_test(y, c("u878", "u867"), gamma = 3)
  expect_equal(ours$deviates[["u878"]], unname(ref$dev["U878"]), tolerance = 0.01)
  expect_equal(ours$p.value, ref$pval, tolerance = 0.01)
})

test_that("adaptive p-value never rejects more often than its level (non-nested regions)", {
  ## 11 untied pairs: the cutoffs move from (5, 10) at 15/256 to (4, 11) at 16/256 so
  ## "smallest rejecting level" alone gives P(p <= 0.0625) = 0.078125
  n <- 11
  sg <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  for (gam in c(1, 1.5)) {
    kap <- gam / (1 + gam)
    w <- apply(sg, 1, function(s) prod(ifelse(s > 0, kap, 1 - kap)))
    ps <- apply(sg, 1, function(s) riposte_adaptive_pvalue(s * seq_len(n), gamma = gam))
    for (a in c(0.02, 0.05, 0.0625, 0.1, 0.2, 0.5))
      expect_lte(sum(w[ps <= a + 1e-12]), a + 1e-12)
  }
})

test_that("adaptive p-value equals the smallest rejecting level when regions are nested", {
  y <- c(5, 4, 3, 2, 1, -0.5, 6, 7, 8, 9, 10, 11)
  p <- riposte_adaptive_pvalue(y)
  expect_true(riposte_adaptive_test(y, alpha = p + 1e-9)$reject)
  expect_false(riposte_adaptive_test(y, alpha = p * 0.99)$reject)
})

test_that("all-zero score columns (six pairs) are dropped, not NaN", {
  set.seed(1); y <- rnorm(6) + 1
  b <- riposte_pairs_test(y, gamma = 1)
  expect_equal(b$p.value, 1); expect_true(all(is.na(b$deviates)))
  H <- cbind(a = rank(abs(y)), z = 0, b = rank(abs(y))^2)
  b2 <- riposte_pairs_bound(y, H, gamma = 1) # one zero column among live ones
  expect_false(is.na(b2$p.value)); expect_true(is.na(b2$deviates[["z"]]))
  pt <- riposte_pairs_perm_test(y, H, nresample = NULL)
  expect_false(is.nan(pt$p.value)); expect_gt(pt$p.value, 0)
  expect_equal(riposte_pairs_perm_test(y, matrix(0, 6, 2), nresample = NULL)$p.value, 1)
})
