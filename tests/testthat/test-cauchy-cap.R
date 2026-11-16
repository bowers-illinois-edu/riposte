## test-cauchy-cap.R
##
## THE POINT. A p-value of exactly 1 converts to minus infinity in Liu and Xie's
## Cauchy combination and sets the combined p-value to 1, whatever the other
## p-values say. riposte's large-sample Cauchy and hybrid avoid that with the
## truncation of Gui, Jiang, and Wang at 0.9 (test-truncated-cauchy.R). The paper
## (Bowers and Burton, footnote in Section 4.5) instead replaces any p-value above
## 0.999 by 0.999 and then applies Liu and Xie's combination. Jake decided on
## 2026-11-16 that riposte offers the paper's rule beside the truncation, which
## stays the default, so that a reader can reproduce the paper's numbers. These
## tests check that the cap is the paper's rule, written out and as the paper's
## own code computes it, and that riposte_test() applies it when asked.

## the paper's rule written out from its definition
capped_by_hand <- function(p, cap = 0.999) {
  stats::pcauchy(mean(tan((0.5 - pmin(p, cap)) * pi)), lower.tail = FALSE)
}

test_that("riposte_capped_cauchy() is the paper's rule", {
  ## the seven p-values in the comment above CAUCHY_P_CAP in the paper's code:
  ## uncapped the combination is 1, capped it is 0.0030
  p <- c(0.0003, 0.15, 0.57, 0.05, 0.81, 0.25, 1)
  expect_equal(riposte_capped_cauchy(p), capped_by_hand(p))
  expect_equal(round(riposte_capped_cauchy(p), 4), 0.0030)
  ## a cap of 1 is Liu and Xie's combination, which a p-value of 1 sets to 1
  expect_equal(riposte_capped_cauchy(p, cap = 1), 1)
  ## below the cap nothing changes
  q <- c(0.01, 0.2, 0.7)
  expect_equal(riposte_capped_cauchy(q), stats::pcauchy(mean(tan((0.5 - q) * pi)),
                                                       lower.tail = FALSE))
  expect_error(riposte_capped_cauchy(p, cap = 0.3), "cap")
})

test_that("riposte_capped_cauchy() agrees with the paper's cauchy_combine()", {
  paper_file <- file.path(Sys.getenv("RIPOSTE_PAPER_REPO", "~/repos/block_test_power"),
                          "Analysis", "simlib", "R", "pvalues.R")
  skip_if_not(file.exists(path.expand(paper_file)), "block_test_power not present")
  paper <- new.env()
  for (f in c("representations", "minp", "pvalues"))
    sys.source(path.expand(file.path(dirname(paper_file), paste0(f, ".R"))), envir = paper)
  set.seed(11)
  for (i in 1:5) {
    p <- c(stats::runif(5), 1, 0.9995)
    expect_equal(riposte_capped_cauchy(p), paper$cauchy_combine(log(p), log1p(-p)))
  }
})

test_that("riposte_test(cauchy_cap = ) applies the paper's rule to the large-sample Cauchy and hybrid", {
  ## Wallsten and Nteta's 36 respondents: the max distance's treated sum sits at
  ## its mean, so its large-sample p-value is exactly 1
  wal <- wallsten_36()
  cau <- riposte_test(Y ~ Z, data = wal, statistic = "cauchy", engine = "asymptotic",
                      cauchy_cap = 0.999)
  expect_equal(cau$p.value, capped_by_hand(cau$component_p))
  hyb <- riposte_test(Y ~ Z, data = wal, statistic = "hybrid", engine = "asymptotic",
                      cauchy_cap = 0.999)
  expect_equal(length(hyb$component_p), 7L)
  expect_equal(hyb$p.value, capped_by_hand(hyb$component_p))
  ## the default is still the truncation at 0.9
  tr <- riposte_test(Y ~ Z, data = wal, statistic = "cauchy", engine = "asymptotic")
  expect_equal(tr$p.value, riposte_truncated_cauchy(tr$component_p, 0.9))
  expect_output(print(cau), "0.999")
})

test_that("riposte_test() takes the cap or the truncation, not both", {
  wal <- wallsten_36()
  expect_error(riposte_test(Y ~ Z, data = wal, statistic = "cauchy", engine = "asymptotic",
                            cauchy_cap = 0.999, cauchy_truncation = 0.8),
               "cauchy_cap")
  expect_error(riposte_test(Y ~ Z, data = wal, statistic = "cauchy", engine = "asymptotic",
                            cauchy_cap = 1.2), "cauchy_cap")
})
