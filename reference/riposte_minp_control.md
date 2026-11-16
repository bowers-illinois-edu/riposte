# Settings for the min-p combination

The min-p combination (`riposte_test(statistic = "minp")`, the default)
takes the smallest of several p-values at the observed assignment and at
every re-randomization, and reports the share of assignments whose
smallest p-value is at most the observed one. This function sets which
p-values it combines and how each is computed.

## Usage

``` r
riposte_minp_control(
  pvalue = c("midp", "randomization", "asymptotic"),
  quadratic = NULL
)
```

## Arguments

- pvalue:

  how each p-value inside the min-p is computed at every assignment.
  `"midp"` (default): the mid-p value of the statistic among the
  assignments, the share more extreme plus half the share tied.
  `"randomization"`: the ordinary randomization p-value, the share of
  assignments at which the statistic is at least as extreme, ties
  counted in full. `"asymptotic"`: the two-sided normal tail of each
  standardized statistic (one-sided under a one-sided `alternative`) and
  the chi-square tail of the quadratic combination, as in the paper's
  code. Under every rule the min-p's own p-value comes from the
  re-randomizations, so the test keeps its level. With `"randomization"`
  the min-p's p-value is also at most K times the smallest of its K
  p-values, so the min-p rejects at level alpha whenever one of them is
  at or below alpha / K. With `"midp"` that bound becomes 2K times the
  smallest mid-p value.

- quadratic:

  whether the quadratic combination's p-value is one of the p-values
  combined. `NULL` (default) includes it for a two-sided test and leaves
  it out for a one-sided test. The quadratic combination squares each
  representation's departure from its re-randomization mean, so it
  responds to departures in both directions; asking for it with a
  one-sided `alternative` is an error.

## Value

a list of class `riposte_minp_control` with elements `pvalue` and
`quadratic`, for the `minp` argument of
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md).

## Examples

``` r
riposte_minp_control()
#> $pvalue
#> [1] "midp"
#> 
#> $quadratic
#> NULL
#> 
#> attr(,"class")
#> [1] "riposte_minp_control"
riposte_minp_control(pvalue = "randomization")
#> $pvalue
#> [1] "randomization"
#> 
#> $quadratic
#> NULL
#> 
#> attr(,"class")
#> [1] "riposte_minp_control"
riposte_minp_control(quadratic = FALSE)
#> $pvalue
#> [1] "midp"
#> 
#> $quadratic
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "riposte_minp_control"
```
