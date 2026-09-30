## test-representations.R
##
## These tests pin down what each representation MEANS, not just that it returns
## a number of the right length. The statistical content: the distance
## representations must equal their textbook definitions, the numerically stable
## form must agree with the naive form to machine precision even when the outcome
## carries a large offset (the improvement over the C++ that motivated rewriting
## these), and the representation set must be open --- a user can supply any
## function of the outcome and have it flow through.

test_that("mean_dist equals the textbook mean L1 distance", {
  set.seed(1)
  y <- rnorm(20)
  ## reference: full distance matrix, row means with denominator n - 1
  ref <- rowSums(abs(outer(y, y, "-"))) / (length(y) - 1)
  expect_equal(riposte_mean_dist(y), ref)
})

test_that("mean_rank_dist is mean_dist applied to mid-ranks", {
  set.seed(2)
  y <- rnorm(15)
  expect_equal(riposte_mean_rank_dist(y),
               riposte_mean_dist(rank(y, ties.method = "average")))
})

test_that("max_dist is the distance to the farther extreme", {
  y <- c(-2, 0, 1, 5)
  d <- riposte_max_dist(y)
  expect_equal(d, pmax(y - min(y), max(y) - y))
  ## each unit's farther extreme is at most the range away; the two extreme units
  ## attain the range (each is a full range from the opposite end)
  expect_equal(max(d), diff(range(y)))
  expect_setequal(which(d == max(d)), c(1L, 4L))
})

test_that("rank is the mid-rank (ties averaged)", {
  y <- c(3, 1, 1, 2)
  expect_equal(riposte_rank(y), rank(y, ties.method = "average"))
})

test_that("tanh representation is bounded and stays informative on a large scale", {
  set.seed(3)
  y <- rnorm(50) * 1000 + 1e6      # an outcome on a large scale (e.g. dollars)
  z <- riposte_tanh(y)
  expect_true(all(abs(z) <= 1))
  ## the point of robust standardisation: it does NOT collapse to +/-1 the way
  ## raw tanh(y) would on this scale (which would carry no information)
  expect_true(stats::sd(z) > 0.1)
  expect_true(all(abs(tanh(y)) > 1 - 1e-8))   # raw tanh is saturated here
})

## --- Huber's psi: the sixth representation of the paper's default set ----------
## The paper (block_test_power, Section "six scores") replaced the tanh transform
## with Huber's psi of the block-standardized outcome: the departure from the
## block median in units of the block's MAD (scaled by 1.4826, as stats::mad()
## does), held within +/- 1.345. The tests below pin down that definition with
## numbers small enough to check by hand, the fallback the paper specifies when
## the MAD is zero, and the property that motivates standardizing first: the
## representation does not depend on the outcome's units or location.

test_that("huber equals the capped, MAD-standardized departure from the median", {
  ## median 3; absolute departures 2, 1, 0, 1, 97, whose median is 1, so the
  ## MAD is 1.4826. Departures in MAD units: -1.349, -0.674, 0, 0.674, 65.4.
  ## Both ends exceed 1.345 in absolute value, so both are capped.
  y <- c(1, 2, 3, 4, 100)
  expect_equal(riposte_huber(y),
               c(-1.345, -1 / 1.4826, 0, 1 / 1.4826, 1.345))
})

test_that("huber caps an outlier rather than letting it dominate, and keeps the order", {
  set.seed(7)
  y <- c(rnorm(29), 50)
  h <- riposte_huber(y)
  expect_true(all(abs(h) <= 1.345))
  expect_equal(h[30], 1.345)                         # capped, not removed
  expect_true(all(diff(h[order(y)]) >= 0))           # nondecreasing in y
})

test_that("huber falls back to the mean absolute deviation when the MAD is zero", {
  ## More than half the block ties at the median 0, so the MAD is 0. The mean
  ## absolute deviation from the median is (5 + 10) / 6 = 2.5, which puts the two
  ## untied units at 2 and 4. With the default cap both are held at 1.345; a cap
  ## of 10 shows the divisor itself.
  y <- c(0, 0, 0, 0, 5, 10)
  expect_equal(riposte_huber(y, k = 10), c(0, 0, 0, 0, 2, 4))
  expect_equal(riposte_huber(y), c(0, 0, 0, 0, 1.345, 1.345))
})

test_that("huber gives zeros to a block with no spread", {
  expect_equal(riposte_huber(rep(2, 6)), rep(0, 6))
  expect_equal(riposte_huber(3), 0)                  # a block of one unit
})

test_that("huber returns an empty vector for an empty block", {
  ## riposte_centred_scores() splits by the block factor, so a block level with
  ## no rows (left over after subsetting a data frame) reaches every
  ## representation as numeric(0). median() of that is NA, so an unguarded
  ## comparison with zero would stop the whole test.
  expect_equal(riposte_huber(numeric(0)), numeric(0))
  y <- c(1, 4, 2, 8, 5, 7)
  block <- factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b", "unused"))
  sm <- riposte_score_matrix(y, block)
  expect_true("huber" %in% sm$kept)
})

test_that("huber does not depend on the outcome's location or units", {
  ## The reason to standardize within block: an outcome in dollars around a
  ## million carries the same information as the same outcome in thousands
  ## around zero. Raw tanh(y) would be saturated at 1 for every unit here.
  set.seed(8)
  y <- rnorm(40)
  expect_equal(riposte_huber(1e6 + 1000 * y), riposte_huber(y))
})

test_that("huber matches the paper's own implementation", {
  ## The paper computes its sixth representation with rep_huber() in
  ## Analysis/simlib/R/representations.R of the block_test_power repository.
  ## That repository is private, so this cross-check runs only where a local
  ## copy exists (set RIPOSTE_PAPER_REPO to point elsewhere).
  paper_file <- file.path(Sys.getenv("RIPOSTE_PAPER_REPO", "~/repos/block_test_power"),
                          "Analysis", "simlib", "R", "representations.R")
  skip_if_not(file.exists(paper_file), "block_test_power not present")
  paper <- new.env()
  sys.source(path.expand(paper_file), envir = paper)
  set.seed(9)
  blocks <- list(normal = rnorm(25),
                 skewed = rexp(18),
                 mad_zero = c(rep(1, 7), 2, 5, 9),
                 constant = rep(4, 5))
  for (nm in names(blocks)) {
    expect_equal(riposte_huber(blocks[[nm]]), paper$rep_huber(blocks[[nm]]),
                 info = nm)
  }
  expect_setequal(names(riposte_reps_default()), paper$SIX)
})

## --- the numerical-stability regression test (the reason these were rewritten) -
## The honest claim: the fast closed form matches the definitional O(n^2)
## sum-of-distances computed on the SAME input to machine precision at any
## location offset. (It cannot recover precision lost when the outcome is stored
## as a large offset plus a small spread --- that floor is in the input doubles,
## not the algorithm --- so we compare to the definition on the same offset x,
## not to the distances of the un-offset data.)
outer_mean_dist <- function(x) rowSums(abs(outer(x, x, "-"))) / (length(x) - 1)

test_that("mean_dist equals the definition on the same input at any offset", {
  set.seed(4)
  base <- rnorm(60)
  for (offset in c(0, 1e6, 1e9, 1e12)) {
    y <- base + offset
    expect_equal(riposte_mean_dist(y), outer_mean_dist(y), tolerance = 1e-9,
                 info = paste("offset =", offset))
  }
})

## manytestsr is an optional, local-only cross-check (it is not on CRAN or any
## remote), so we reach its C++ function by namespace rather than declaring a
## formal dependency; both tests skip when manytestsr is not installed (e.g. CI).
mt_mean_dist <- function(y) {
  get("fast_dists_and_trans_hybrid", asNamespace("manytestsr"))(y)$mean_dist
}

test_that("mean_dist adds less cancellation than the prefix-sum C++ at a large offset", {
  skip_if_not_installed("manytestsr")
  set.seed(5)
  base <- rnorm(60)
  y <- base + 1e9                       # a realistic large location
  def <- outer_mean_dist(y)             # the definition on the same input
  err_riposte <- max(abs(riposte_mean_dist(y) - def))
  err_cpp <- max(abs(mt_mean_dist(y) - def))
  ## our centred form should be orders of magnitude closer to the definition
  expect_lt(err_riposte, 1e-9)
  expect_gt(err_cpp, 1e-9)
  expect_lt(err_riposte, err_cpp)
})

test_that("mean_dist matches manytestsr's value when there is no offset", {
  skip_if_not_installed("manytestsr")
  set.seed(6)
  y <- rnorm(30)
  expect_equal(riposte_mean_dist(y), mt_mean_dist(y), tolerance = 1e-10)
})

## --- extensibility: representations are an open, user-supplied set -------------
test_that("the default set is the paper's six, with Huber's psi as the sixth", {
  reps <- riposte_reps_default()
  expect_setequal(names(reps),
                  c("raw", "rank", "mean_dist", "mean_rank_dist", "max_dist", "huber"))
  expect_true(all(vapply(reps, is.function, logical(1))))
  expect_identical(reps$huber, riposte_huber)
})

test_that("the proposal set keeps six representations, with raw tanh as the sixth", {
  ## The NSF proposal (and manytestsr's pIndepDist) used tanh of the raw outcome
  ## in the sixth place. Swapping it in must replace Huber's psi, not add a
  ## seventh representation.
  reps <- riposte_reps_proposal()
  expect_setequal(names(reps),
                  c("raw", "rank", "mean_dist", "mean_rank_dist", "max_dist", "tanh"))
  y <- c(-2000, 0, 0.5, 3000)
  expect_equal(reps$tanh(y), tanh(y))
})

test_that("a user-supplied transformation flows through the score matrix", {
  set.seed(6)
  y <- rnorm(12)
  block <- factor(rep(1:2, each = 6))
  ## a custom representation: signed square root, applied within block
  my_reps <- c(riposte_reps_default(),
               list(signed_sqrt = function(v) sign(v) * sqrt(abs(v))))
  sm <- riposte_score_matrix(y, block, representations = my_reps)
  expect_true("signed_sqrt" %in% sm$kept)
  ## the kept column equals the within-block-centred custom score
  raw <- sign(y) * sqrt(abs(y))
  centred <- raw - ave(raw, block)
  expect_equal(sm$scores[, "signed_sqrt"], centred, ignore_attr = TRUE)
})

test_that("validate_reps rejects malformed representation lists", {
  expect_error(riposte_validate_reps(list()), "non-empty")
  expect_error(riposte_validate_reps(list(function(y) y)), "unique")          # unnamed
  expect_error(riposte_validate_reps(list(a = 1)), "function")
})
