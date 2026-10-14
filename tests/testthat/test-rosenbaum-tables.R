## Full replication of the numeric tables in Rosenbaum (2012), Biometrika 99:763-774
## (Tables 2 and 3) from the 679 smoker/non-smoker pairs. The data ship with the CRAN
## package `testtwice` (data set `smokerlead`), so these tests are skipped without it.
## The paper prints three decimals, so we allow a tolerance of 0.0011 (it also appears
## to truncate in a few cells; see the notes below).

smoker_scores <- function() {
  skip_if_not_installed("testtwice"); skip_if_not_installed("mvtnorm")
  e <- new.env()
  utils::data("smokerlead", package = "testtwice", envir = e)
  d <- e$smokerlead
  stopifnot(all(sign(d$lead) == sign(d$llead)))     # the paper combines lead and log-lead
  sc <- function(y, type, ...) riposte_pairs_scores(y, type, ...)
  list(y = d$lead, dose = d$dose,
       W_nd = sc(d$lead, "wilcoxon"), W_d = sc(d$lead, "wilcoxon") * d$dose,
       B_nd = sc(d$lead, "brown"),    B_d = sc(d$lead, "brown") * d$dose,
       U878 = sc(d$lead, "ustat", 8, 7, 8), U867 = sc(d$lead, "ustat", 8, 6, 7),
       lW_nd = sc(d$llead, "wilcoxon"), lW_d = sc(d$llead, "wilcoxon") * d$dose,
       lB_nd = sc(d$llead, "brown"),    lB_d = sc(d$llead, "brown") * d$dose,
       lU878 = sc(d$llead, "ustat", 8, 7, 8), lU867 = sc(d$llead, "ustat", 8, 6, 7))
}

## upper bound on the one-sided p-value for ONE statistic (normal approximation)
one_p <- function(y, h, gamma) {
  th <- gamma / (1 + gamma)
  1 - stats::pnorm((sum((y > 0) * h) - th * sum(h)) / sqrt(th * (1 - th) * sum(h^2)))
}

test_that("Biometrika Table 2: 12 single-statistic sensitivity analyses", {
  s <- smoker_scores()
  cols_lead <- c("W_nd", "W_d", "B_nd", "B_d", "U878", "U867")
  cols_log  <- c("lW_nd", "lW_d", "lB_nd", "lB_d", "lU878", "lU867")
  lead <- rbind(  # Gamma = 1, 2, 2.3, 2.5, 2.8, 3 ; columns as cols_lead
    c(0, 0, 0, 0, 0, 0), c(0, .001, 0, .001, 0, 0), c(.002, .035, 0, .020, .006, 0),
    c(.016, .131, .004, .075, .026, .001), c(.147, .429, .045, .273, .119, .015),
    c(.351, .652, .133, .462, .236, .050))
  g_lead <- c(1, 2, 2.3, 2.5, 2.8, 3)
  logs <- rbind(  # Gamma = 1, 2, 2.5, 2.8, 3, 3.2, 3.4, 3.6 ; columns as cols_log
    c(0, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 0, 0), c(0, .023, 0, .005, 0, 0),
    c(.012, .139, .001, .036, 0, .001), c(.052, .299, .004, .093, .001, .004),
    c(.154, .496, .016, .192, .003, .014), c(.324, .684, .048, .326, .009, .041),
    c(.531, .826, .113, .478, .022, .095))
  g_log <- c(1, 2, 2.5, 2.8, 3, 3.2, 3.4, 3.6)
  for (i in seq_along(g_lead)) for (j in seq_along(cols_lead))
    expect_lte(abs(one_p(s$y, s[[cols_lead[j]]], g_lead[i]) - lead[i, j]), 0.0011,
               label = paste("lead", cols_lead[j], "Gamma", g_lead[i]))
  for (i in seq_along(g_log)) for (j in seq_along(cols_log))
    expect_lte(abs(one_p(s$y, s[[cols_log[j]]], g_log[i]) - logs[i, j]), 0.0011,
               label = paste("log", cols_log[j], "Gamma", g_log[i]))
})

test_that("Biometrika Table 3: joint procedure (6) for L = 3, 3, 4, 10, 12", {
  s <- smoker_scores()
  sets <- list(
    Wilcoxon = c("W_nd", "W_d", "lW_nd"), Brown = c("B_nd", "B_d", "lB_nd"),
    Ustat = c("U878", "U867", "lU878", "lU867"),
    All10 = c("W_nd", "W_d", "lW_nd", "B_nd", "B_d", "lB_nd", "U878", "U867", "lU878", "lU867"),
    All12 = c("W_nd", "W_d", "B_nd", "B_d", "U878", "U867",
              "lW_nd", "lW_d", "lB_nd", "lB_d", "lU878", "lU867"))
  tab <- rbind(  # Gamma = 1, 2, 2.5, 2.8, 3, 3.2, 3.4, 3.6
    c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0), c(.001, 0, 0, 0, 0), c(.020, .001, 0, .001, .001),
    c(.084, .007, .003, .004, .004), c(.224, .029, .008, .012, .013),
    c(.431, .082, .022, .030, .033), c(.647, .180, .049, .065, .069))
  gam <- c(1, 2, 2.5, 2.8, 3, 3.2, 3.4, 3.6)
  set.seed(20261013)       # L = 10 and 12 use mvtnorm's Monte Carlo algorithm
  for (i in seq_along(gam)) for (j in seq_along(sets)) {
    H <- do.call(cbind, s[sets[[j]]])
    p <- riposte_pairs_bound(s$y, H, gamma = gam[i])$p.value
    ## Notes: two cells of the paper differ from ours by less than 0.001 and look like
    ## truncation: Brown at Gamma = 3 (ours 0.0076, paper .007) and U-statistic at
    ## Gamma = 2.8 (ours 0.0006, paper .000).
    expect_lte(abs(p - tab[i, j]), 0.0011,
               label = paste(names(sets)[j], "Gamma", gam[i]))
  }
})

test_that("Annals Proposition 1: design sensitivities for Normal effects, tau/sigma = 1/2", {
  ## Computed from the paper's Proposition 1 for Normal effects with tau/sigma = 1/2;
  ## the paper reports 4.97 (Noether) and 3.60 (Brown). zeta(lambda) is the chance an
  ## effect-shifted draw falls in the top lambda band of |y|.
  zeta <- function(lam, tau = 0.5) {
    f <- function(cc) stats::pnorm(cc - tau) - stats::pnorm(-cc - tau) - (1 - lam)
    cc <- stats::uniroot(f, c(0, 20))$root
    1 - stats::pnorm(cc - tau)
  }
  z1 <- zeta(1 / 3); z2 <- zeta(2 / 3)
  expect_equal(round(z1 / (1 / 3 - z1), 2), 4.97)
  expect_equal(round((z1 + z2) / ((1 / 3 - z1) + (2 / 3 - z2)), 2), 3.60)
})
