## spike-acat-size.R --- does the ANALYTIC (Liu-Xie) ACAT combined p-value
## control the false-positive rate in riposte's unadjusted block-randomized
## setting, and if not, what drives the inflation?
##
## THE QUESTION (Jake). riposte calibrates the ACAT statistic by PERMUTATION, not
## by Liu-Xie's analytic standard-Cauchy tail. The recollection is that the
## analytic version had false-positive-rate (FPR) problems. We want to (a) SEE
## whether that is true in the unadjusted case, and (b) LEARN the mechanism, so
## the package design rests on evidence, not memory. Covariance adjustment is OUT
## OF SCOPE here (a deeper structural problem, set aside on purpose).
##
## METHOD: a diagnostic ladder that turns the candidate mechanisms on one at a
## time, so any inflation can be attributed rather than just observed.
##   Rung 1  independent, EXACTLY uniform marginals  -> is the analytic tail itself
##           calibrated? (Liu-Xie proved yes here; validates our implementation.)
##   Rung 2  DEPENDENT, exactly uniform marginals     -> does dependence alone break
##           it? compound-symmetry sweep + riposte's actual representation
##           correlation. Isolates DEPENDENCE.
##   Rung 3  full riposte: dependent AND discrete permutation marginals, analytic
##           vs the permutation calibration. The added gap over Rung 2 is the
##           DISCRETENESS/non-uniformity of the marginals; the permutation column
##           is the size-valid reference.
## Liu & Xie (2020, JASA 115:393-402): the ACAT statistic sum_j w_j tan((0.5-p_j)pi)
## (weights summing to 1) has an approximately standard-Cauchy tail under H0, so
## the analytic p-value is 0.5 - atan(T)/pi.

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(20260625L)

## analytic Liu-Xie combined p-value from a matrix of marginal p-values
## (rows = datasets, cols = representations); equal weights w_j = 1/d
analytic_acat_p <- function(P) {
  term <- tan((0.5 - P) * pi)
  sm <- P < 1e-15; term[sm] <- 1 / (P[sm] * pi)            # pole guards (rarely fire)
  bg <- P > 1 - 1e-15; term[bg] <- -1 / ((1 - P[bg]) * pi)
  0.5 - atan(rowMeans(term)) / pi
}
fpr    <- function(p, a) mean(p <= a)
fpr_se <- function(a, n) sqrt(a * (1 - a) / n)              # MC SE under the nominal rate
ALPHA  <- c(0.10, 0.05, 0.01)

## two-sided p-values that are EXACTLY uniform under H0, from latent normals Z;
## dependence enters only through cor(Z) -- a Gaussian copula with uniform margins
unif_marginals_from_Z <- function(Z) 2 * stats::pnorm(-abs(Z))

cat("=== Does analytic ACAT control the FPR (unadjusted)? ===\n")
cat(sprintf("(MC SE at alpha=.05 is +/- %.4f for Nsim=5e4, +/- %.4f for Nsim=2500)\n",
            fpr_se(.05, 5e4), fpr_se(.05, 2500)))

## ---- Rung 1: independent, exactly uniform marginals --------------------------
Nbig <- 50000L; d <- 6L
P1 <- matrix(runif(Nbig * d), Nbig, d)
p_an1 <- analytic_acat_p(P1)
cat("\n[Rung 1] independent, exactly-uniform marginals (analytic ACAT's ideal):\n")
for (a in ALPHA) cat(sprintf("  alpha=%.2f  FPR=%.4f\n", a, fpr(p_an1, a)))

## ---- Rung 2: dependence only (uniform marginals, Gaussian copula) -------------
## (a) compound-symmetry sweep: one shared factor drives correlation rho
cat("\n[Rung 2a] dependent uniform marginals, compound symmetry (isolates dependence):\n")
for (rho in c(0.0, 0.3, 0.6, 0.9)) {
  Fac <- matrix(rnorm(Nbig), Nbig, d)                       # idiosyncratic
  shared <- rnorm(Nbig)
  Z <- sqrt(rho) * shared + sqrt(1 - rho) * Fac             # cor(Z_j, Z_k) = rho
  p_an <- analytic_acat_p(unif_marginals_from_Z(Z))
  cat(sprintf("  rho=%.1f  FPR(.10)=%.4f  FPR(.05)=%.4f  FPR(.01)=%.4f\n",
              rho, fpr(p_an, .10), fpr(p_an, .05), fpr(p_an, .01)))
}

## (b) riposte's ACTUAL representation correlation, with uniform marginals.
## Build it from one representative null dataset's Strasser-Weber Sigma.
mk_design <- function(n_blocks, bs, m, effect = 0) {
  n <- n_blocks * bs
  block <- factor(rep(seq_len(n_blocks), each = bs))
  z <- numeric(n); for (ix in split(seq_len(n), block)) z[ix[seq_len(m)]] <- 1
  y <- rnorm(n) + effect * z
  list(y = y, z = z, block = block, n = n)
}
rep_dn <- mk_design(40L, 6L, 3L)
Sref <- riposte_score_matrix(rep_dn$y, rep_dn$block)
Mref <- riposte_sw_moments(Sref$scores, rep_dn$z, rep_dn$block)
Rrep <- stats::cov2cor(Mref$Sigma)
cat("\n[Rung 2b] dependence = riposte's actual representation correlation, uniform marginals:\n")
cat("  representations:", paste(Sref$kept, collapse = ", "), "\n")
cat(sprintf("  mean |off-diagonal correlation| = %.2f, max = %.2f\n",
            mean(abs(Rrep[lower.tri(Rrep)])), max(abs(Rrep[lower.tri(Rrep)]))))
L <- chol(Rrep)                                             # Z = E %*% L has cor = Rrep
dref <- ncol(Rrep)
Z <- matrix(rnorm(Nbig * dref), Nbig, dref) %*% L
p_an2b <- analytic_acat_p(unif_marginals_from_Z(Z))
for (a in ALPHA) cat(sprintf("  alpha=%.2f  FPR=%.4f\n", a, fpr(p_an2b, a)))

## ---- Rung 3: full riposte, analytic vs permutation, across block size --------
## Now the marginals are riposte's two-sided mid-p PERMUTATION p-values: dependent
## AND discrete. Fix y, redraw the observed assignment each iteration (the actual
## randomization under H0), compute both combined p-values on identical data.
## Smaller blocks -> more discrete marginals, so the trend across bs isolates the
## discreteness contribution; the permutation column should hold size throughout.
cat("\n[Rung 3] full riposte: analytic vs PERMUTATION calibration, by block size:\n")
cat("  (marginals = two-sided mid-p; permutation column is the size-valid reference)\n")
Nsim <- 2500L; B <- 999L
run_riposte_fpr <- function(bs, m, n_blocks = 40L) {
  dn <- mk_design(n_blocks, bs, m)                          # fixed y, block
  scores <- riposte_score_matrix(dn$y, dn$block)$scores
  idx <- split(seq_len(dn$n), dn$block)
  p_an <- p_pm <- numeric(Nsim)
  for (i in seq_len(Nsim)) {
    z <- numeric(dn$n)                                      # one fresh assignment
    for (ix in idx) z[ix[sample.int(length(ix), m)]] <- 1
    res <- riposte_cauchy(scores, z, dn$block, nresample = B)
    p_pm[i] <- res$p.value                                  # permutation of T_c
    p_an[i] <- analytic_acat_p(matrix(res$component_midp, 1))
  }
  list(p_an = p_an, p_pm = p_pm)
}
designs <- list(c(bs = 4L, m = 2L), c(bs = 6L, m = 3L), c(bs = 12L, m = 6L))
cat(sprintf("  %-22s %18s %18s\n", "design", "analytic FPR", "permutation FPR"))
for (de in designs) {
  r <- run_riposte_fpr(de["bs"], de["m"])
  cat(sprintf("  bs=%2d m=%d (C=%2d/blk)   .10/.05/.01: %.3f/%.3f/%.3f   %.3f/%.3f/%.3f\n",
              de["bs"], de["m"], choose(de["bs"], de["m"]),
              fpr(r$p_an, .10), fpr(r$p_an, .05), fpr(r$p_an, .01),
              fpr(r$p_pm, .10), fpr(r$p_pm, .05), fpr(r$p_pm, .01)))
}
cat(sprintf("\n  (Nsim=%d, B=%d; MC SE at .05 ~ %.4f)\n", Nsim, B, fpr_se(.05, Nsim)))
