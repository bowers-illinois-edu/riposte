# Choose (or blend) the combination from the conditioning of Sigma

Given a permutation covariance `Sigma` of the representations, returns
the screen's decision. Under `"threshold"`, if the condition number of
the correlation matrix exceeds `threshold`, choose the Cauchy (which
inverts nothing), else the quadratic. Under `"shrink"` (the default),
shrink the correlation matrix toward the identity,
`R_a = (1 - a) R + a I`, choosing the intensity `a` so the condition
number drops to `threshold`; `a = 0` is the full quadratic and `a = 1`
is the diagonal covariance (which ignores the cross-representation
dependence, as the Cauchy does). The shrunk covariance is returned for
the quadratic to use.

## Usage

``` r
riposte_screen(Sigma, control = riposte_screen_control())
```

## Arguments

- Sigma:

  a permutation covariance of the representations (the closed-form
  covariance without adjustment, the pooled empirical covariance with
  it).

- control:

  a
  [`riposte_screen_control()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_screen_control.md)
  list.

## Value

a list with `method`, the `condition` number of the correlation matrix,
and either `choice` (`"quadratic"`/`"cauchy"`, for `"threshold"`) or the
shrinkage `lambda`, the shrunk `Sigma`, and the achieved
`condition_shrunk` (for `"shrink"`).

## Details

Choosing the combination from the condition number does not disturb the
test's exact level, but the reason differs by path. Without covariance
adjustment, `Sigma` is the closed-form permutation covariance — a
function of the outcomes (fixed under the sharp null) and the design,
NOT of the realized assignment — so the condition number is constant
across the randomization reference, an ancillary quantity. With
covariance adjustment, `Sigma` is the pooled empirical covariance of the
residual statistics, which depends on the realized draws; it is no
longer ancillary, but it is a symmetric function of the exchangeable
pooled set, which is the property exactness actually needs (see
`dev/theory-adjusted-exactness.md`).

The shrinkage-intensity rule here is a sensible default, not the settled
Project 1 result; expect it to change. Set `lambda` in
[`riposte_screen_control()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_screen_control.md)
to fix the intensity yourself. Exact level holds for any shrinkage rule;
the rule affects power only.
