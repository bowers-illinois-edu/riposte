## acat-helpers.R --- shared helpers for the ACAT false-positive-rate spikes.
## Sourced by spike-acat-size.R and spike-acat-size-stress.R. No package code.

## analytic Liu-Xie combined p-value from a matrix of marginal p-values
## (rows = datasets, cols = representations); equal weights w_j = 1/d.
## Liu & Xie (2020): sum_j w_j tan((0.5-p_j)pi) has an approx standard-Cauchy
## tail under H0, so the combined p-value is 0.5 - atan(T)/pi.
analytic_acat_p <- function(P) {
  term <- tan((0.5 - P) * pi)
  sm <- P < 1e-15; term[sm] <- 1 / (P[sm] * pi)            # pole guards (rarely fire)
  bg <- P > 1 - 1e-15; term[bg] <- -1 / ((1 - P[bg]) * pi)
  0.5 - atan(rowMeans(term)) / pi
}

fpr    <- function(p, a) mean(p <= a)
fpr_se <- function(a, n) sqrt(a * (1 - a) / n)              # MC SE under the nominal rate

## two-sided p-values that are EXACTLY uniform under H0 from latent normals Z;
## dependence enters only through cor(Z) (a Gaussian copula with uniform margins).
## Note: tan((0.5 - p)pi) of such a p is exactly standard Cauchy, so each ACAT
## term is marginally standard Cauchy and only the JOINT law carries dependence.
unif_marginals_from_Z <- function(Z) 2 * stats::pnorm(-abs(Z))

## one block-randomized design: fixed y (optional sharp effect), m treated/block
mk_design <- function(n_blocks, bs, m, effect = 0) {
  n <- n_blocks * bs
  block <- factor(rep(seq_len(n_blocks), each = bs))
  z <- numeric(n); for (ix in split(seq_len(n), block)) z[ix[seq_len(m)]] <- 1
  y <- stats::rnorm(n) + effect * z
  list(y = y, z = z, block = block, n = n)
}
