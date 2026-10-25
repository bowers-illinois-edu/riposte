## level-truncated-cauchy.R --- rejection rates of the large-sample truncated
## Cauchy and hybrid combinations under the sharp null, with many more draws
## than the level tests in tests/testthat can afford. Run from the package
## root: Rscript dev/level-truncated-cauchy.R. Results of the 2026-10-25 run:
##
##   Wallsten and Nteta's 36 outcomes held fixed, 100,000 reassignments
##   (SE 0.0007 at 0.05, 0.0003 at 0.01):
##              cauchy t=0.9  cauchy t=1  hybrid t=0.9  hybrid t=1
##     at 0.05     0.0487       0.0436       0.0471       0.0415
##     at 0.01     0.0063       0.0058       0.0050       0.0046
##
##   Fresh 5-point outcomes, 36 units, half treated, 50,000 data sets
##   (SE 0.0010 at 0.05, 0.0004 at 0.01):
##     at 0.05     0.0507       0.0467       0.0495       0.0454
##     at 0.01     0.0077       0.0070       0.0073       0.0066
##
## For contrast, with 40 independent p-values the truncated combination
## rejects at about 0.063 at 0.05 (RItools, vignettes/truncated-cauchy-memo.md).

## --- 1. Wallsten and Nteta's 36 outcomes, treatment reassigned
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-wallsten.R")
d0 <- wallsten_36()
RNGkind("L'Ecuyer-CMRG"); set.seed(20261025)
B <- 100000; cores <- 8
chunks <- split(seq_len(B), rep(seq_len(cores), length.out = B))
res <- parallel::mclapply(chunks, function(idx) {
  t(vapply(idx, function(i) {
    d <- d0; d$Z <- sample(d$Z)
    c(cauchy_t09 = riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic")$p.value,
      cauchy_t1  = riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic", cauchy_truncation = 1)$p.value,
      hybrid_t09 = riposte_test(Y ~ Z, d, statistic = "hybrid", engine = "asymptotic")$p.value,
      hybrid_t1  = riposte_test(Y ~ Z, d, statistic = "hybrid", engine = "asymptotic", cauchy_truncation = 1)$p.value)
  }, numeric(4)))
}, mc.cores = cores, mc.set.seed = TRUE)
P <- do.call(rbind, res)
saveRDS(P, "dev/level-truncated-cauchy-wallsten.rds")
rates <- rbind(at05 = colMeans(P <= 0.05), at01 = colMeans(P <= 0.01))
print(round(rates, 4))
cat("B =", nrow(P), " SE at 0.05:", round(sqrt(.05*.95/nrow(P)), 5), " SE at 0.01:", round(sqrt(.01*.99/nrow(P)), 5), "\n")

## --- 2. fresh data sets with no effect

RNGkind("L'Ecuyer-CMRG"); set.seed(20261026)
B <- 50000; cores <- 12
chunks <- split(seq_len(B), rep(seq_len(cores), length.out = B))
res <- parallel::mclapply(chunks, function(idx) {
  t(vapply(idx, function(i) {
    d <- data.frame(Y = sample(c(0, 0.25, 0.5, 0.75, 1), 36, replace = TRUE), Z = sample(rep(0:1, 18)))
    c(cauchy_t09 = riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic")$p.value,
      cauchy_t1  = riposte_test(Y ~ Z, d, statistic = "cauchy", engine = "asymptotic", cauchy_truncation = 1)$p.value,
      hybrid_t09 = riposte_test(Y ~ Z, d, statistic = "hybrid", engine = "asymptotic")$p.value,
      hybrid_t1  = riposte_test(Y ~ Z, d, statistic = "hybrid", engine = "asymptotic", cauchy_truncation = 1)$p.value)
  }, numeric(4)))
}, mc.cores = cores, mc.set.seed = TRUE)
P2 <- do.call(rbind, res)
saveRDS(P2, "dev/level-truncated-cauchy-random.rds")
rates <- rbind(at05 = colMeans(P2 <= 0.05), at01 = colMeans(P2 <= 0.01))
print(round(rates, 4))
cat("B =", nrow(P2), " SE at 0.05:", round(sqrt(.05*.95/nrow(P2)), 5), " SE at 0.01:", round(sqrt(.01*.99/nrow(P2)), 5), "\n")
