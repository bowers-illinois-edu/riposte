## test-clusters.R
##
## Cluster randomization. The statistical content:
##   - the collapse is correct (cluster means; treatment and block carried;
##     validation of the two cluster-randomization facts);
##   - the cluster test equals a block-randomized test run on the collapsed data
##     (the reduction the implementation rests on);
##   - CLUSTERING MATTERS: when treatment is cluster-assigned and outcomes are
##     correlated within cluster, analyzing at the unit level (ignoring clusters)
##     is badly anti-conservative, while the cluster test holds the level. This is
##     the reason the extension exists.
##   - exactness at the cluster level (full enumeration of cluster assignments).

make_cluster_data <- function(seed = 1, effect = 0, B = 3L, G = 4L, m = 4L) {
  set.seed(seed)
  ncl <- B * G
  block_of_cluster <- rep(seq_len(B), each = G)
  clust <- rep(seq_len(ncl), each = m)
  block <- rep(block_of_cluster, each = m)
  u <- rnorm(ncl, sd = 1.5)                       # cluster random effect
  zc <- integer(ncl)
  for (b in seq_len(B)) { ix <- which(block_of_cluster == b); zc[ix][sample.int(G, G %/% 2L)] <- 1L }
  z <- zc[clust]
  y <- u[clust] + rnorm(length(clust), sd = 0.7) + effect * z
  data.frame(Y = y, trt = z, blk = factor(block), cl = clust)
}

test_that("collapse produces cluster means and validates the cluster facts", {
  d <- make_cluster_data()
  col <- riposte:::riposte_collapse_clusters(d$Y, d$trt, d$blk, d$cl)
  ## cluster-level outcome is the cluster mean
  expect_equal(col$y, as.numeric(tapply(d$Y, d$cl, mean)), ignore_attr = TRUE)
  expect_equal(col$n_clusters, length(unique(d$cl)))
  expect_equal(col$n_units, nrow(d))
  ## treatment varying within a cluster is rejected
  bad <- d; bad$trt[1] <- 1 - bad$trt[1]
  expect_error(riposte:::riposte_collapse_clusters(bad$Y, bad$trt, bad$blk, bad$cl),
               "treatment varies within cluster")
  ## a cluster spanning two blocks is rejected
  bad2 <- d; bad2$blk[1] <- levels(bad2$blk)[length(levels(bad2$blk))]
  expect_error(riposte:::riposte_collapse_clusters(bad2$Y, bad2$trt, bad2$blk, bad2$cl),
               "must nest within blocks")
})

test_that("the cluster test equals a block test on the collapsed data", {
  d <- make_cluster_data(effect = 1)
  ## manual collapse to one row per cluster
  agg <- function(v, by) as.numeric(tapply(v, by, function(x) x[1]))
  cd <- data.frame(
    Y = as.numeric(tapply(d$Y, d$cl, mean)),
    trt = agg(d$trt, d$cl),
    blk = factor(agg(as.integer(d$blk), d$cl))
  )
  a <- riposte_test(Y ~ trt | blk, d, clusters = "cl", statistic = "max",
                    nresample = 199, seed = 7)
  b <- riposte_test(Y ~ trt | blk, cd, statistic = "max", nresample = 199, seed = 7)
  expect_equal(a$p.value, b$p.value)
  expect_equal(a$statistic, b$statistic)
})

test_that("the cluster test reports clusters and the effective n", {
  d <- make_cluster_data()
  r <- riposte_test(Y ~ trt | blk, d, clusters = "cl", statistic = "max",
                    nresample = 99, seed = 1)
  expect_true(isTRUE(r$clustered))
  expect_equal(r$n_clusters, length(unique(d$cl)))
  expect_equal(r$n, length(unique(d$cl)))           # effective n = number of clusters
  expect_output(print(r), "clusters")
})

test_that("the cluster test is exact under full enumeration of cluster assignments", {
  ## 1 block, 4 clusters of 2 units, 2 clusters treated -> choose(4,2)=6 assignments
  set.seed(3)
  clust <- rep(1:4, each = 2)
  y <- rnorm(8)
  col <- riposte:::riposte_collapse_clusters(y, rep(c(1, 1, 0, 0), each = 2),
                                             factor(rep(1, 8)), clust)
  block_c <- col$block
  G <- enumerate_assignments(block_c)               # cluster-level enumeration (6 cols)
  sm <- riposte_score_matrix(col$y, block_c)
  mom <- riposte_sw_moments(sm$scores, G[, 1], block_c)
  Tmat <- crossprod(G, sm$scores); cen <- sweep(Tmat, 2, mom$mu)
  q <- rowSums((cen %*% MASS::ginv(mom$Sigma)) * cen)
  p_exact <- mean(q >= q[1])
  p_rip <- riposte_quadratic(sm$scores, G[, 1], block_c, moments = mom, draws = G)$p.value
  expect_equal(p_rip, p_exact)
})

test_that("ignoring clustering is anti-conservative; the cluster test holds the level", {
  skip_on_cran()
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260623)
  ## Calibrated (n = 120, ICC high): naive unit-level size ~0.57, cluster ~0.067.
  B <- 4L; G <- 6L; m <- 5L; ncl <- B * G
  block_of_cluster <- rep(seq_len(B), each = G)
  clust <- rep(seq_len(ncl), each = m); block <- rep(block_of_cluster, each = m)
  size <- function(nsims, clustered) {
    p <- numeric(nsims)
    for (i in seq_len(nsims)) {
      u <- rnorm(ncl, sd = 2)
      y <- u[clust] + rnorm(length(clust), sd = 0.7)   # sharp null
      zc <- integer(ncl)
      for (b in seq_len(B)) { ix <- which(block_of_cluster == b); zc[ix][sample.int(G, G %/% 2L)] <- 1L }
      d <- data.frame(Y = y, trt = zc[clust], blk = factor(block), cl = clust)
      p[i] <- if (clustered)
        riposte_test(Y ~ trt | blk, d, clusters = "cl", statistic = "max", nresample = 99)$p.value
      else
        riposte_test(Y ~ trt | blk, d, statistic = "max", nresample = 99)$p.value
    }
    mean(p <= 0.05)
  }
  expect_gt(size(100L, FALSE), 0.20)                  # naive: badly inflated
  cluster_size <- size(100L, TRUE)
  expect_lt(cluster_size, 0.12)                       # cluster test: holds the level
  expect_gt(cluster_size, 0.01)
})

test_that("covariance adjustment works under cluster randomization", {
  d <- make_cluster_data(effect = 1)
  set.seed(2); d$x1 <- rnorm(nrow(d)); d$x2 <- rnorm(nrow(d))
  r <- riposte_test(Y ~ trt | blk, d, clusters = "cl", statistic = "max",
                    representations = list(raw = function(y) y),
                    adjust = ~ x1 + x2, nresample = 99, seed = 1)
  expect_true(isTRUE(r$clustered) && isTRUE(r$adjusted))
  expect_true(r$p.value > 0 && r$p.value <= 1)
})
