## test-minp-block.R
##
## THE POINT. Every p-value in the min-p over seven compares treated with
## control units within each block and adds the comparisons over blocks, so a
## treatment that helps in some blocks and harms in others can cancel to
## nothing. The block statistic squares each block's standardized comparison
## before adding, so no block's sign cancels another's. Bowers and Burton add
## it, with the cosine score, as an eighth p-value for block-randomized
## experiments. Jake decided on 2026-11-16 (PLAN_block_statistic_2026-11-16.md)
## that riposte offers it for any score through riposte_minp_control(block = ),
## off by default until the paper chooses a score. These tests check:
##   1. the statistic is its formula, with the chi-square on the number of
##      blocks in the sum as its large-sample reference;
##   2. what the formula implies at the two extremes: with one block it is the
##      square of the score's standardized statistic, and with blocks of two
##      every block's term is 1 at every assignment;
##   3. inside the min-p, under every rule, riposte equals a brute-force count,
##      the test is exact, and with randomization p-values the bound K times
##      the smallest still holds with K = 8;
##   4. with the cosine score and the asymptotic rule riposte gives the paper's
##      own numbers;
##   5. it finds help in half the blocks and harm in the rest, which the rank
##      alone misses;
##   6. riposte_test() offers it only where it can add something, and says so.

rules <- c("midp", "randomization", "asymptotic")

## three blocks of four units, two treated in each: 6^3 = 216 assignments
three_block_enumeration <- function(seed = 20261117) {
  set.seed(seed)
  block <- factor(rep(1:3, each = 4))
  y <- rnorm(12)
  G <- enumerate_assignments(block)
  sm <- riposte_score_matrix(y, block)
  Tmat <- crossprod(G, sm$scores)
  colnames(Tmat) <- colnames(sm$scores)
  list(block = block, y = y, G = G, scores = sm$scores, Tmat = Tmat,
       rank = within_block_rank(y, block))
}

## ---- the statistic -------------------------------------------------------------

test_that("the block statistic is its formula on a full enumeration", {
  e <- three_block_enumeration()
  bf <- brute_block_stat(e$rank, e$block, e$G)
  res <- riposte:::riposte_block_squared(e$rank, e$block, e$G)
  expect_equal(res$stat, bf$stat, tolerance = 1e-10)
  expect_equal(res$df, 3L)
  ## its null mean over every assignment is the number of blocks: each term is a
  ## squared standardized sum, whose mean is its variance divided by itself
  expect_equal(mean(res$stat), 3, tolerance = 1e-10)
})

test_that("with one block the block statistic is the square of the standardized statistic", {
  set.seed(5)
  block <- factor(rep(1, 10))
  y <- rnorm(10)
  G <- enumerate_assignments(block)
  score <- rank(y)
  res <- riposte:::riposte_block_squared(score, block, G)
  mom <- enum_moments(crossprod(G, matrix(score - mean(score))))
  t_all <- drop(crossprod(G, score - mean(score)))
  expect_equal(res$stat, (t_all - mom$mu)^2 / drop(mom$Sigma), tolerance = 1e-10)
  expect_equal(res$df, 1L)
})

test_that("with blocks of two every block's term is 1 at every assignment", {
  ## in a block of two with one treated, the treated departure from the block
  ## mean is plus or minus half the difference, and its variance is the square
  ## of that half, so the term is 1 whichever unit is treated
  set.seed(6)
  block <- factor(rep(1:5, each = 2))
  y <- rnorm(10)
  G <- enumerate_assignments(block)
  res <- riposte:::riposte_block_squared(y, block, G)
  expect_equal(res$stat, rep(5, ncol(G)), tolerance = 1e-10)
})

test_that("a block whose scores all tie, or that has one arm only, is left out", {
  block <- factor(rep(1:3, each = 4))
  score <- c(1, 2, 3, 4,   5, 5, 5, 5,   1, 3, 2, 9)   # block 2 ties throughout
  z <- c(1, 1, 0, 0,   1, 0, 1, 0,   1, 1, 0, 0)
  G <- riposte_block_draws(z, block, 49L)
  res <- riposte:::riposte_block_squared(score, block, G)
  expect_equal(res$df, 2L)
  expect_equal(res$stat, brute_block_stat(score, block, G)$stat, tolerance = 1e-10)
})

## ---- inside the min-p --------------------------------------------------------

test_that("the min-p over eight p-values equals a brute-force count under every rule", {
  e <- three_block_enumeration()
  bs <- brute_block_stat(e$rank, e$block, e$G)
  for (rule in rules) {
    P <- cbind(brute_member_p(e$Tmat, rule, quadratic = TRUE),
               block_rank = brute_block_p(bs, rule))
    expect_equal(ncol(P), 8L)
    p <- brute_minp_all(P)
    res <- riposte:::riposte_minp(e$scores, e$G[, 1], e$block, draws = e$G,
                                  pvalue = rule, quadratic = TRUE,
                                  block_scores = cbind(rank = e$rank))
    expect_equal(res$p.value, p[1], label = rule)
    expect_equal(res$component_p, P[1, ], tolerance = 1e-8, label = rule)
  }
})

test_that("the min-p over eight p-values is exact, and the K-times bound holds with K = 8", {
  e <- three_block_enumeration()
  bs <- brute_block_stat(e$rank, e$block, e$G)
  for (rule in rules) {
    P <- cbind(brute_member_p(e$Tmat, rule, quadratic = TRUE),
               block_rank = brute_block_p(bs, rule))
    p <- brute_minp_all(P)
    levels <- sort(unique(p))
    excess <- vapply(levels, function(a) mean(p <= a + 1e-12) - a, numeric(1))
    expect_true(all(excess <= 1e-12), label = paste(rule, "level"))
    if (rule == "randomization")
      expect_true(all(p <= 8 * apply(P, 1, min) + 1e-12))
  }
})

## ---- the paper's own numbers --------------------------------------------------

test_that("with the cosine score and the asymptotic rule riposte gives the paper's min-p over eight", {
  ## The paper computes the cosine statistic with member_square_then_sum() and
  ## the min-p over eight with p_minp_eight() in Analysis/simlib/R/minp.R of the
  ## block_test_power repository, which is private; this runs only where a
  ## local copy exists (set RIPOSTE_PAPER_REPO to point elsewhere).
  paper_dir <- file.path(Sys.getenv("RIPOSTE_PAPER_REPO", "~/repos/block_test_power"),
                         "Analysis", "simlib", "R")
  skip_if_not(dir.exists(path.expand(paper_dir)), "block_test_power not present")
  paper <- new.env()
  for (f in c("representations", "minp", "pvalues"))
    sys.source(path.expand(file.path(paper_dir, paste0(f, ".R"))), envir = paper)
  ## base-R replacement for the paper's data.table step, as in test-minp-rules.R
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

  d <- help_and_harm_design(20261117, B = 6L, nb = 10L)
  block <- factor(d$b)
  ## riposte's cosine score is the paper's
  expect_equal(unlist(lapply(split(d$y, block), riposte_cosine)),
               unlist(lapply(split(d$y, block), function(v) paper$score_cos(v, 1))))

  out <- list(dt = list(Y = d$y, blockF = block), Z = d$z)
  RNGkind("L'Ecuyer-CMRG"); set.seed(8)
  G <- riposte_block_draws(d$z, block, 999L)
  P_paper <- cbind(paper$member_singles(out)(G), paper$member_quadratic(out)(G),
                   paper$member_square_then_sum(out)(G))
  smallest <- apply(P_paper, 1, min)
  p_paper <- (sum(smallest[-1] <= smallest[1] + 1e-12) + 1) / ncol(G)

  sm <- riposte_score_matrix(d$y, block)
  cosine <- stats::ave(d$y, block, FUN = riposte_cosine)
  res <- riposte:::riposte_minp(sm$scores, d$z, block, draws = G, pvalue = "asymptotic",
                                quadratic = TRUE, block_scores = cbind(cosine = cosine))
  expect_equal(res$p.value, p_paper)
  expect_equal(res$component_p[["block_cosine"]], P_paper[1, 8], tolerance = 1e-8)
})

## ---- what it finds ---------------------------------------------------------------

test_that("the block statistic finds help in half the blocks and harm in the rest", {
  ## ten blocks of twenty; treatment raises outcomes by 0.8 in five blocks and
  ## lowers them by 0.8 in the other five, so the rank's comparisons cancel
  d <- help_and_harm_design(1)
  res <- riposte_test(y ~ z | b, d, nresample = 999, seed = 1,
                      minp = riposte_minp_control(block = "rank"))
  expect_gt(res$component_p[["rank"]], 0.5)
  expect_lt(res$component_p[["block_rank"]], 0.01)
  expect_lt(res$p.value, 0.05)
})

## ---- riposte_test() ------------------------------------------------------------

test_that("riposte_test() adds a block statistic only when asked", {
  d <- help_and_harm_design(2, B = 4L, nb = 8L)
  none <- riposte_test(y ~ z | b, d, nresample = 199, seed = 2)
  expect_equal(length(none$component_p), 7L)
  expect_null(none$minp$block)

  by_name <- riposte_test(y ~ z | b, d, nresample = 199, seed = 2,
                          minp = riposte_minp_control(block = "rank"))
  expect_equal(names(by_name$component_p)[8], "block_rank")
  expect_equal(by_name$minp$block, "rank")

  by_function <- riposte_test(y ~ z | b, d, nresample = 199, seed = 2,
                              minp = riposte_minp_control(block = list(cosine = riposte_cosine)))
  expect_equal(names(by_function$component_p)[8], "block_cosine")
  expect_output(print(by_function), "block statistic")
})

test_that("riposte_minp_control() and riposte_test() check the block setting", {
  expect_error(riposte_minp_control(block = 3), "block")
  expect_error(riposte_minp_control(block = list(function(y) y)), "named")
  d <- help_and_harm_design(3, B = 4L, nb = 8L)
  expect_error(riposte_test(y ~ z | b, d, nresample = 49,
                            minp = riposte_minp_control(block = "cosine")),
               "representations")
})

test_that("riposte_test() leaves the block statistic out where it cannot add anything", {
  ## one block: the square of the score's own standardized statistic
  d1 <- help_and_harm_design(4, B = 1L, nb = 20L)
  expect_message(
    one <- riposte_test(y ~ z, d1, nresample = 99, seed = 4,
                        minp = riposte_minp_control(block = "rank")),
    "left out")
  expect_equal(length(one$component_p), 7L)
  ## blocks of two: every block's term is 1, so the statistic never changes
  d2 <- help_and_harm_design(5, B = 10L, nb = 2L)
  expect_message(
    two <- riposte_test(y ~ z | b, d2, nresample = 99, seed = 5,
                        minp = riposte_minp_control(block = "rank")),
    "left out")
  expect_false("block_rank" %in% names(two$component_p))
})

test_that("a one-sided or covariance-adjusted min-p with a block statistic stops with a reason", {
  d <- help_and_harm_design(6, B = 4L, nb = 8L)
  msg <- tryCatch(
    riposte_test(y ~ z | b, d, alternative = "greater", nresample = 49,
                 minp = riposte_minp_control(block = "rank")),
    error = conditionMessage)
  expect_match(msg, "block statistic", fixed = TRUE)
  expect_match(msg, "both directions", fixed = TRUE)
  d$x <- rnorm(nrow(d))
  expect_error(riposte_test(y ~ z | b, d, adjust = ~ x, nresample = 49,
                            minp = riposte_minp_control(block = "rank")),
               "covariance adjustment")
})

test_that("riposte_cosine() is the cosine of the within-block rank scaled into (0, 1)", {
  y <- c(3.2, -1, 0.5, 7, 0.5)
  r <- rank(y)
  expect_equal(riposte_cosine(y), cos(pi * r / (length(y) + 1)))
  ## near 1 at the lowest outcome and near -1 at the highest
  expect_gt(riposte_cosine(y)[which.min(y)], 0.8)
  expect_lt(riposte_cosine(y)[which.max(y)], -0.8)
})
