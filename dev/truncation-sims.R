## truncation-sims.R --- the simulations behind vignettes/truncated-cauchy.Rmd.
##
## THE POINT. riposte's large-sample Cauchy combination converts each p-value
## p to tan((0.5 - t p) pi) and refers the sum S of the converted values to
## min(1, n P(Cauchy > S) / t), following Gui, Jiang, and Wang (2025). The
## vignette explains why t = 0.9. These simulations measure, for several values
## of t, the false positive rate when the treatment changes nothing and the
## power when it shifts or spreads the outcome, and they compare one value of t
## with the smallest p-value over several. Results are summarized in
## inst/extdata/truncation_sims.rds, which the vignette reads, because the
## simulations take several minutes.
##
## Run from the package root:  Rscript dev/truncation-sims.R
## Adapted from the scripts of two reviews written for the Sarkar and Coppock
## reanalysis in the block_test_power repository (Analysis/truncation_review/).
suppressMessages(devtools::load_all(".", quiet = TRUE))
RNGkind("L'Ecuyer-CMRG"); set.seed(20261005)

## The 36 respondents of Wallsten and Nteta (2016), comparison 6, as counts of
## each answer by arm (Sarkar and Coppock 2026 replication archive). The
## max distance's treated mean equals its overall mean, 31/36, so its
## large-sample p-value is exactly 1.
vals <- c(0, 0.25, 0.5, 0.75, 1)
wal <- data.frame(Y = c(rep(vals, c(10, 2, 3, 2, 1)), rep(vals, c(5, 2, 0, 8, 3))),
                  Z = rep(0:1, each = 18))

## The six large-sample component p-values for many assignments at once. The
## scores depend only on Y, and the exact mean and covariance of the treated
## sums only on the design, so one matrix product gives every assignment's
## sums. Checked against riposte_test() below.
comp_many <- function(y, Zmat) {
  X <- riposte_score_matrix(y, rep(1, length(y)))$scores
  mom <- riposte_sw_moments(X, Zmat[, 1], rep(1, length(y)))
  cen <- sweep(crossprod(Zmat, X), 2, mom$mu)
  return(pchisq(sweep(cen^2, 2, diag(mom$Sigma), "/"), df = 1, lower.tail = FALSE))
}
z_chk <- sample(rep(0:1, 20)); y_chk <- rlnorm(40)
stopifnot(isTRUE(all.equal(
  as.numeric(comp_many(y_chk, cbind(z_chk))),
  as.numeric(riposte_test(Y ~ Z, data = data.frame(Y = y_chk, Z = z_chk), statistic = "cauchy",
                          engine = "asymptotic")$component_p))))

## Combining rules applied row by row to a matrix of component p-values.
trunc_p <- function(Pm, t) pmin(1, ncol(Pm) * pcauchy(rowSums(1 / tan(pi * t * Pm)), lower.tail = FALSE) / t)
rules <- list(
  `untruncated (Liu and Xie)` = function(Pm) apply(Pm, 1, riposte_truncated_cauchy, truncation = 1),
  `t = 0.95` = function(Pm) trunc_p(Pm, 0.95), `t = 0.9` = function(Pm) trunc_p(Pm, 0.9),
  `t = 0.8` = function(Pm) trunc_p(Pm, 0.8), `t = 0.7` = function(Pm) trunc_p(Pm, 0.7),
  `t = 0.5` = function(Pm) trunc_p(Pm, 0.5),
  `harmonic mean (t near 0)` = function(Pm) pmin(1, ncol(Pm) / rowSums(1 / Pm)))
## the trunc_p shortcut must agree with the package function
stopifnot(isTRUE(all.equal(trunc_p(comp_many(y_chk, cbind(z_chk)), 0.9),
                           riposte_truncated_cauchy(as.numeric(comp_many(y_chk, cbind(z_chk))), 0.9))))

## 1. False positive rate and power for each rule, all rules applied to the
## same component p-values. False positives: outcomes held fixed, treatment
## reassigned 20,000 times. Power: 4,000 experiments with 40 lognormal
## outcomes, half treated, the treatment adding 0.55 or multiplying treated
## outcomes' distance from the control median by 2.1.
R0 <- 20000L; R1 <- 4000L; z40 <- rep(0:1, 20)
y_ln <- rlnorm(40)
P <- list(
  `false positives, Wallsten` = comp_many(wal$Y, replicate(R0, sample(wal$Z))),
  `false positives, lognormal` = comp_many(y_ln, replicate(R0, sample(z40))),
  `power, shift` = t(replicate(R1, { y0 <- rlnorm(40); z <- sample(z40); as.numeric(comp_many(y0 + 0.55 * z, cbind(z))) })),
  `power, spread` = t(replicate(R1, { y0 <- rlnorm(40); z <- sample(z40); m <- median(y0[z == 0])
    as.numeric(comp_many(ifelse(z == 1, m + 2.1 * (y0 - m), y0), cbind(z))) })))
rates <- sapply(P, function(Pm) sapply(rules, function(r) mean(r(Pm) <= 0.05)))
draws <- sapply(P, nrow)

## 2. Six independent uniform p-values, the case Gui, Jiang, and Wang's
## Figure 2 finds hardest for these rules at 0.05.
U <- matrix(runif(6 * 1e6), ncol = 6)
independent <- sapply(rules[-1], function(r) mean(r(U) <= 0.05))
independent <- c(`untruncated (Liu and Xie)` = mean(pcauchy(rowMeans(tan((0.5 - U) * pi)), lower.tail = FALSE) <= 0.05),
                 independent)

## 3. One truncation level against the smallest p-value over several, with
## p-values from the large-sample formula and from re-randomization. For
## re-randomization, each experiment's statistic (its combined p-value) is
## compared with the same statistic on 999 reassignments.
B <- 999L; N <- 2000L
several <- function(y, z) {
  Pm <- comp_many(y, cbind(z, replicate(B, sample(z))))
  stat <- cbind(`t = 0.9` = trunc_p(Pm, 0.9), `smallest over t = 0.5, 0.7, 0.9` =
                  pmin(trunc_p(Pm, 0.5), trunc_p(Pm, 0.7), trunc_p(Pm, 0.9)))
  perm <- apply(stat, 2, function(s) (1 + sum(s[-1] <= s[1])) / (B + 1))
  return(c(large_sample = stat[1, ] <= 0.05, rerandomization = perm <= 0.05))
}
multi <- list(
  `false positives, Wallsten` = t(replicate(N, several(wal$Y, sample(wal$Z)))),
  `power, shift` = t(replicate(N, { y0 <- rlnorm(40); z <- sample(z40); several(y0 + 0.55 * z, z) })),
  `power, spread` = t(replicate(N, { y0 <- rlnorm(40); z <- sample(z40); m <- median(y0[z == 0])
    several(ifelse(z == 1, m + 2.1 * (y0 - m), y0), z) })))
several_rates <- sapply(multi, colMeans)

out <- list(rates = rates, draws = draws, independent = independent, several = several_rates,
            several_draws = N, B = B, date = Sys.Date(), riposte = as.character(packageVersion("riposte")))
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(out, "inst/extdata/truncation_sims.rds")
print(out)
