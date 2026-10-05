## helper-wallsten.R --- a small coarse data set from the Sarkar and Coppock
## reanalysis, used by test-hybrid.R and test-truncated-cauchy.R.

## The 36 respondents of Wallsten and Nteta (2016), comparison 6, by arm. The
## treated respondents' mean of max(y, 1 - y) is 31/36, exactly the mean over
## all 36, so the max distance's treated sum sits at its mean.
wallsten_36 <- function() {
  vals <- c(0, 0.25, 0.5, 0.75, 1)
  data.frame(Y = c(rep(vals, c(10, 2, 3, 2, 1)), rep(vals, c(5, 2, 0, 8, 3))),
             Z = rep(0:1, each = 18))
}
