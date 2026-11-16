## test-batched-draws.R
##
## THE POINT. Every re-randomization statistic in riposte needs, at each redrawn
## assignment, only one number per representation (and one per block statistic).
## Until 0.0.0.9015 riposte held every redrawn assignment in memory at once: an
## n by (nresample + 1) matrix of doubles, about 32 GB at 100,000 units and
## 40,000 redraws. Jake's task of 2026-11-06 asked riposte to draw in batches.
## riposte now draws a batch of assignments, keeps each one's statistics, and
## discards the batch, so memory grows with the batch rather than with
## nresample. The option riposte.batch_entries sets the largest batch, in matrix
## entries (default 1e7, about 80 MB). These tests check:
##   1. the batches are the assignments riposte_block_draws() makes, in the same
##      order from the same seed, the observed assignment first, and no batch is
##      larger than the option allows;
##   2. riposte_test() gives the same result whatever the batch size, for every
##      statistic, with a block statistic, with covariance adjustment, and with
##      clusters, so batching changes memory and nothing else;
##   3. an experiment of 100,000 units runs.

batch_design <- function(seed, B = 4L, nb = 10L) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  block <- factor(rep(seq_len(B), each = nb))
  z <- integer(B * nb)
  for (b in levels(block)) {
    ix <- which(block == b); z[ix][sample.int(nb, nb %/% 2L)] <- 1L
  }
  data.frame(y = rnorm(B * nb) + 0.5 * z, z = z, b = block, x = rnorm(B * nb))
}

## the same kind of experiment with treatment assigned to clusters of two units
clustered_batch_design <- function(seed, B = 4L, clusters_per_block = 6L) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  cl_block <- rep(seq_len(B), each = clusters_per_block)
  cl_z <- unlist(lapply(seq_len(B), function(b)
    sample(rep(0:1, each = clusters_per_block / 2))))
  unit_cl <- rep(seq_along(cl_block), each = 2L)
  data.frame(y = rnorm(length(unit_cl)) + 0.5 * cl_z[unit_cl], z = cl_z[unit_cl],
             b = factor(cl_block[unit_cl]), cl = factor(unit_cl))
}

test_that("the batches are riposte_block_draws()'s assignments, in order, none too large", {
  d <- batch_design(1)
  block <- factor(d$b)
  RNGkind("L'Ecuyer-CMRG"); set.seed(9)
  G <- riposte_block_draws(d$z, block, 57L)
  widths <- integer(0)
  RNGkind("L'Ecuyer-CMRG"); set.seed(9)
  cols <- riposte:::riposte_over_draws(d$z, block, 57L, function(Gb) {
    widths <<- c(widths, ncol(Gb))
    t(Gb)
  }, entries = 10 * nrow(d))
  expect_equal(t(cols), G, ignore_attr = TRUE)
  expect_equal(sum(widths), 58L)
  expect_true(all(widths <= 10L))
  expect_equal(cols[1, ], as.numeric(d$z))
})

test_that("riposte_test() gives the same result whatever the batch size", {
  d <- batch_design(2)
  run <- function(...) riposte_test(y ~ z | b, d, nresample = 199, seed = 2, ...)
  calls <- list(
    minp = list(),
    minp_block = list(minp = riposte_minp_control(block = "rank")),
    minp_rand = list(minp = riposte_minp_control(pvalue = "randomization")),
    quadratic = list(statistic = "quadratic"),
    cauchy = list(statistic = "cauchy"),
    max = list(statistic = "max"),
    hybrid = list(statistic = "hybrid"),
    screen = list(statistic = "screen"),
    one_sided = list(alternative = "greater"),
    adjusted = list(adjust = ~ x),
    adjusted_fit_once = list(adjust = ~ x, adjust_refit = FALSE)
  )
  dc <- clustered_batch_design(2)
  run_cl <- function() riposte_test(y ~ z | b, dc, clusters = "cl", nresample = 199, seed = 2)
  whole <- lapply(calls, function(a) do.call(run, a))
  whole$clustered <- run_cl()
  old <- options(riposte.batch_entries = 7 * nrow(d))
  on.exit(options(old))
  pieces <- lapply(calls, function(a) do.call(run, a))
  pieces$clustered <- run_cl()
  calls$clustered <- list()
  for (nm in names(calls)) {
    expect_equal(pieces[[nm]]$p.value, whole[[nm]]$p.value, label = nm)
    expect_equal(pieces[[nm]]$statistic, whole[[nm]]$statistic, label = nm)
    expect_equal(pieces[[nm]]$component_p, whole[[nm]]$component_p, label = nm)
  }
})

test_that("riposte_components() gives the same result whatever the batch size", {
  d <- batch_design(3)
  whole <- riposte_components(y ~ z | b, d, nresample = 199, seed = 3)
  old <- options(riposte.batch_entries = 5 * nrow(d))
  on.exit(options(old))
  pieces <- riposte_components(y ~ z | b, d, nresample = 199, seed = 3)
  expect_equal(pieces$component_p, whole$component_p)
})

test_that("an experiment of 100,000 units runs", {
  skip_on_cran()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(4)
  n <- 100000L
  block <- factor(rep(seq_len(100L), each = 1000L))
  z <- as.integer(stats::ave(seq_len(n), block, FUN = function(i) sample(rep(0:1, length(i) / 2))))
  d <- data.frame(y = rnorm(n) + 0.02 * z, z = z, b = block)
  res <- riposte_test(y ~ z | b, d, nresample = 199, seed = 4)
  expect_true(res$p.value > 0 && res$p.value <= 1)
  expect_equal(res$nresample, 199L)
})
