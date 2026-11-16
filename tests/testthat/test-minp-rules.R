## test-minp-rules.R
##
## THE POINT. The min-p takes the smallest of K p-values at every assignment
## and refers it to that smallest value at the other assignments. Jake decided
## on 2026-11-16 (MEMO_paper_changes_2026-11-16.md) that riposte lets the user
## choose how the K p-values are computed, and that the quadratic combination
## can be one of them, as it is in the paper's min-p over seven p-values.
## These tests pin down what each choice promises:
##   1. the level: under every rule, on a full enumeration, the min-p's p-value
##      falls at or below any level a on at most a share a of the assignments;
##   2. the arithmetic: under every rule, riposte's p-value and its K p-values
##      equal a brute-force count (helper-minp.R);
##   3. the guarantee the paper's next version leads with: with randomization
##      p-values inside, the min-p's p-value is at most K times the smallest of
##      the K, so the min-p rejects at 0.05 whenever one of them is at or below
##      0.05/K. With mid-p values inside, that bound can fail; the weaker bound
##      2K times the smallest mid-p value holds instead;
##   4. the "asymptotic" rule is the paper's own min-p: riposte gives the
##      number that Analysis/simlib/R/minp.R gives on the same assignments.
## Each enumeration uses two blocks of six units, three treated in each (400
## assignments), unless it says otherwise.

rules <- c("midp", "randomization", "asymptotic")

two_block_enumeration <- function(seed = 20261116) {
  set.seed(seed)
  block <- factor(rep(1:2, each = 6))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  sm <- riposte_score_matrix(y, block)
  Tmat <- crossprod(G, sm$scores)
  colnames(Tmat) <- colnames(sm$scores)
  list(block = block, G = G, scores = sm$scores, Tmat = Tmat)
}

## ---- the level, under every rule ---------------------------------------------

test_that("under every rule the min-p over seven p-values is an exact test", {
  ## The smallest of the seven p-values is a function of the assignment once the
  ## outcomes are fixed, and the min-p refers it to the same function at every
  ## assignment the design allows. So with each assignment taken in turn as the
  ## observed one, the min-p's p-value is at most a on at most a share a of
  ## them, whatever rule produced the seven.
  e <- two_block_enumeration()
  for (rule in rules) {
    P <- brute_member_p(e$Tmat, rule, quadratic = TRUE)
    expect_equal(ncol(P), 7L)
    p <- brute_minp_all(P)
    levels <- sort(unique(p))
    excess <- vapply(levels, function(a) mean(p <= a + 1e-12) - a, numeric(1))
    expect_true(all(excess <= 1e-12),
                label = sprintf("%s: share at or below %.4f exceeds it", rule,
                                levels[which.max(excess)]))
  }
})

## ---- the arithmetic, under every rule ----------------------------------------

test_that("riposte's min-p and its seven p-values equal a brute-force count under every rule", {
  e <- two_block_enumeration()
  for (rule in rules) {
    P <- brute_member_p(e$Tmat, rule, quadratic = TRUE)
    p <- brute_minp_all(P)
    ## the observed assignment, and the assignment with the smallest p-value
    for (j in unique(c(1L, which.min(p)))) {
      res <- riposte:::riposte_minp(e$scores, e$G[, j], e$block,
                                    draws = observed_first(e$G, j),
                                    pvalue = rule, quadratic = TRUE)
      expect_equal(res$p.value, p[j], label = paste(rule, "p-value at", j))
      expect_equal(res$component_p, P[j, ], tolerance = 1e-8,
                   label = paste(rule, "seven p-values at", j))
      expect_equal(res$statistic, min(P[j, ]), tolerance = 1e-8)
      expect_true(res$quadratic)
      expect_equal(res$pvalue, rule)
    }
  }
})

test_that("without the quadratic combination the default rule is the six-p-value min-p of before", {
  ## quadratic = FALSE and pvalue = "midp" must reproduce what riposte_minp()
  ## computed before 2026-11-16, which test-minp.R checks against its own count
  e <- two_block_enumeration()
  old <- riposte:::riposte_minp(e$scores, e$G[, 1], e$block, draws = e$G)
  new <- riposte:::riposte_minp(e$scores, e$G[, 1], e$block, draws = e$G,
                                pvalue = "midp", quadratic = FALSE)
  expect_equal(new$p.value, old$p.value)
  expect_equal(new$component_p, old$component_p)
  expect_equal(length(new$component_p), 6L)
  expect_false("quadratic" %in% names(new$component_p))
})

test_that("one-sided rules equal a brute-force count", {
  ## one-sided min-p combines the six one-sided p-values; the quadratic
  ## combination has no one-sided form and stays out
  e <- two_block_enumeration()
  for (alt in c("greater", "less")) for (rule in rules) {
    P <- brute_member_p(e$Tmat, rule, alternative = alt)
    p <- brute_minp_all(P)
    res <- riposte:::riposte_minp(e$scores, e$G[, 1], e$block, draws = e$G,
                                  alternative = alt, pvalue = rule, quadratic = FALSE)
    expect_equal(res$p.value, p[1], label = paste(alt, rule))
    expect_equal(res$component_p, P[1, ], tolerance = 1e-8, label = paste(alt, rule))
  }
})

test_that("every rule measures departures from the null mean it is given", {
  ## Without covariance adjustment the treated sums of within-block-centered
  ## scores have null mean zero. With adjustment riposte refits the residuals
  ## at every assignment, and the null mean is the pooled mean over the
  ## observed and redrawn assignments, which is not zero; the max and the
  ## quadratic combination already measure departures from it. Every rule of
  ## the min-p does the same: shifting every treated sum by a constant, and
  ## the mean with it, leaves every p-value unchanged.
  e <- two_block_enumeration()
  mom <- riposte_sw_moments(e$scores, e$G[, 1], e$block)
  offset <- seq_len(ncol(e$Tmat)) * 0.37
  shifted <- sweep(e$Tmat, 2, offset, "+")
  for (alt in c("two.sided", "greater")) for (rule in rules) {
    quad <- alt == "two.sided"
    a <- riposte:::riposte_minp_from_T(e$Tmat, alt, pvalue = rule, quadratic = quad,
                                       mu = mom$mu, Sigma = mom$Sigma)
    b <- riposte:::riposte_minp_from_T(shifted, alt, pvalue = rule, quadratic = quad,
                                       mu = mom$mu + offset, Sigma = mom$Sigma)
    expect_equal(b$p.value, a$p.value, label = paste(alt, rule))
    expect_equal(b$component_p, a$component_p, tolerance = 1e-8, label = paste(alt, rule))
  }
})

## ---- the guarantee -------------------------------------------------------------

test_that("with randomization p-values the min-p's p-value is at most K times the smallest", {
  ## The proof (Section 4 of the paper's memo of 11 November): an assignment
  ## counts toward the min-p's p-value only if one of its K p-values is at or
  ## below m, the observed smallest; for each test the share of such
  ## assignments is at most m when the p-values are randomization p-values;
  ## adding the K shares gives at most K m.
  e <- two_block_enumeration()
  for (alt in c("two.sided", "greater", "less")) {
    quad <- alt == "two.sided"
    P <- brute_member_p(e$Tmat, "randomization", alternative = alt, quadratic = quad)
    K <- ncol(P)
    p <- brute_minp_all(P)
    expect_true(all(p <= K * apply(P, 1, min) + 1e-12), label = alt)
    ## riposte at the assignment where the bound is tightest
    j <- which.max(p / apply(P, 1, min))
    res <- riposte:::riposte_minp(e$scores, e$G[, j], e$block,
                                  draws = observed_first(e$G, j), alternative = alt,
                                  pvalue = "randomization", quadratic = quad)
    expect_lte(res$p.value, K * min(res$component_p) + 1e-12)
  }
})

test_that("with mid-p values the K-times bound can fail, and the 2K-times bound holds", {
  ## One block of 12 units, 6 treated: 924 assignments. Two tests: the raw
  ## outcome, and the indicator of being above the third quartile, whose treated
  ## sum takes few values. A mid-p value is smaller than the ordinary p-value by
  ## half the share tied with it, so the share of assignments with a mid-p value
  ## at or below m can exceed m. It cannot exceed 2m: where the mid-p value
  ## P - a/2 is at most m, P is at most 2m - (P - a), so at most 2m.
  block <- factor(rep(1, 12))
  G <- enumerate_assignments(block)
  set.seed(1)
  y <- rnorm(12)
  sm <- riposte_score_matrix(y, block, list(raw = function(v) v, top = top_quartile))
  Tmat <- crossprod(G, sm$scores)
  K <- 2L
  P_rand <- brute_member_p(Tmat, "randomization")
  P_mid <- brute_member_p(Tmat, "midp")
  p_rand <- brute_minp_all(P_rand)
  p_mid <- brute_minp_all(P_mid)
  smallest_ordinary <- apply(P_rand, 1, min)

  expect_true(all(p_rand <= K * smallest_ordinary + 1e-12))
  expect_true(all(p_mid <= 2 * K * apply(P_mid, 1, min) + 1e-12))
  ## the mid-p min-p exceeds K times the smallest ordinary p-value somewhere
  ratio <- p_mid / (K * smallest_ordinary)
  expect_gt(max(ratio), 1)

  ## riposte reproduces both at that assignment
  j <- which.max(ratio)
  Gj <- observed_first(G, j)
  mid <- riposte:::riposte_minp(sm$scores, G[, j], block, draws = Gj, pvalue = "midp")
  rand <- riposte:::riposte_minp(sm$scores, G[, j], block, draws = Gj,
                                 pvalue = "randomization")
  expect_equal(mid$p.value, p_mid[j])
  expect_equal(rand$p.value, p_rand[j])
  expect_gt(mid$p.value, K * min(rand$component_p))
  expect_lte(rand$p.value, K * min(rand$component_p) + 1e-12)
})

## ---- the paper's own min-p ------------------------------------------------------

test_that("the asymptotic rule gives the paper's min-p over seven p-values", {
  ## The paper computes its min-p with p_minp_seven() in
  ## Analysis/simlib/R/minp.R of the block_test_power repository: the six
  ## single-representation p-values from the chi-square on one degree of
  ## freedom and the quadratic combination's from the chi-square on the rank
  ## of the covariance, the smallest referred to redrawn assignments. That
  ## repository is private, so this cross-check runs only where a local copy
  ## exists (set RIPOSTE_PAPER_REPO to point elsewhere).
  paper_dir <- file.path(Sys.getenv("RIPOSTE_PAPER_REPO", "~/repos/block_test_power"),
                         "Analysis", "simlib", "R")
  skip_if_not(dir.exists(path.expand(paper_dir)), "block_test_power not present")
  paper <- new.env()
  for (f in c("representations", "minp", "pvalues"))
    sys.source(path.expand(file.path(paper_dir, paste0(f, ".R"))), envir = paper)
  ## The paper's representation_matrix() builds its table with data.table,
  ## which riposte does not use. This replacement does what it does in base R:
  ## apply each of the paper's six representations within blocks, drop one that
  ## is constant within every block, and divide each column by its standard
  ## deviation.
  paper$representation_matrix <- function(out, reps = paper$SIX) {
    Y <- out$dt$Y
    b <- out$dt$blockF
    D <- vapply(reps, function(nm) {
      v <- numeric(length(Y))
      for (ix in split(seq_along(Y), b)) v[ix] <- paper$REPRESENTATIONS[[nm]](Y[ix])
      v
    }, numeric(length(Y)))
    constant <- apply(D, 2, function(v) all(tapply(v, b, function(u) length(unique(u)) < 2L)))
    D <- D[, !constant, drop = FALSE]
    sweep(D, 2, apply(D, 2, stats::sd), "/")
  }

  d <- minp_interface_design(20261116, B = 4L, nb = 10L, effect = 0.5)
  block <- factor(d$b)
  out <- list(dt = list(Y = d$y, blockF = block), Z = d$z)
  RNGkind("L'Ecuyer-CMRG"); set.seed(7)
  G <- riposte_block_draws(d$z, block, 999L)
  P_paper <- cbind(paper$member_singles(out)(G), paper$member_quadratic(out)(G))
  smallest <- apply(P_paper, 1, min)
  p_paper <- (sum(smallest[-1] <= smallest[1] + 1e-12) + 1) / ncol(G)

  sm <- riposte_score_matrix(d$y, block)
  res <- riposte:::riposte_minp(sm$scores, d$z, block, draws = G,
                                pvalue = "asymptotic", quadratic = TRUE)
  expect_equal(res$p.value, p_paper)
  ## the seven observed p-values, matched by name (the paper orders the six
  ## representations differently)
  paper_obs <- stats::setNames(P_paper[1, ], c(paper$SIX, "quadratic"))
  expect_equal(res$component_p[names(paper_obs)], paper_obs, tolerance = 1e-8)
})
