## Fixed outcomes avoid Monte Carlo noise in comparisons of the two programs.
## Block sizes and treated shares differ so that ignoring the design changes
## the null variances. No observed component statistic is exactly zero here.
asymptotic_example <- function() {
  data.frame(
    Y = c(-3, -1, 0, 0.5, 2, 4, 7, 9,
          -2, -1, 0, 1, 2, 3, 4, 5, 8, 13,
          -5, -2, -1, 0, 1, 1.5, 3, 6, 8, 9, 12, 15),
    trt = c(1, 0, 1, 0, 0, 1, 0, 0,
            1, 0, 1, 1, 0, 1, 0, 1, 1, 0,
            1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0),
    blk = factor(rep(1:3, c(8, 10, 12)))
  )
}

## Use coin independently of riposte's score centering and moment calculation.
## This is the interface the paper uses for its quadratic statistic and its
## standardized single-representation statistics.
asymptotic_coin_reference <- function(d, reps = riposte_reps_default()) {
  values <- lapply(reps, function(f) {
    ans <- numeric(nrow(d))
    for (ix in split(seq_len(nrow(d)), d$blk)) ans[ix] <- f(d$Y[ix])
    return(ans)
  })
  dat <- as.data.frame(values)
  dat$trt <- factor(d$trt)
  dat$blk <- d$blk
  fmla <- as.formula(paste(paste(names(reps), collapse = " + "),
                           "~ trt | blk"))
  return(coin::independence_test(fmla, data = dat,
                               teststat = "quadratic",
                               distribution = "asymptotic"))
}
