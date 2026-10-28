# P-value of the exact adaptive test

Write `alpha*` for the smallest level at which the adaptive test rejects
the observed counts. Rosenbaum chooses the cutoffs `(kB, kT)` separately
at each level, so the rejection regions at different levels are not
always nested: counts rejected at one level can go unrejected at a
larger one, and `alpha*` itself, reported as a p-value, would fall at or
below a level more often than that level allows. The reported p-value
treats `alpha*` as a test statistic. It is the probability, under the
sensitivity bound, of counts whose smallest rejecting level is at most
the observed `alpha*`, which is the probability of the union of the
rejection regions at every level up to `alpha*`. It is exact at
`gamma = 1` and an upper bound for `gamma > 1`, so
`P(p.value <= a) <= a`. When the regions are nested it is at most
`alpha*`. Rosenbaum (2012, p. 100) reports `alpha*` itself as the
p-value;
[`riposte_adaptive_smallest_level()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_adaptive_smallest_level.md)
returns it.

## Usage

``` r
riposte_adaptive_pvalue(y, gamma = 1, lambda = c(1/3, 2/3), scan = NULL)
```

## Arguments

- y:

  matched-pair differences (treated minus control). Zero differences are
  dropped. Ties among `|y|` are broken by order
  (`ties.method = "first"`) with a warning, because the paper assumes
  untied ranks.

- gamma:

  sensitivity parameter.

- lambda:

  the two band fractions `0 < lambda[1] < lambda[2] <= 1` (default
  `c(1/3, 2/3)`).

- scan:

  not used; kept so that calls written for an earlier version still run.

## Value

a single p-value in (0, 1\].

## Details

The chance of each rejection region is computed as 1 minus the chance of
the counts outside it, and that subtraction loses all its digits once
the region's chance is near `1e-16`. Below about `1e-15`, rounding error
rather than Rosenbaum's conditions decides the cutoffs. A reported
p-value that small is still valid, because it is the exact chance of the
region the code chose, but it should be read only as "below about
`1e-15`".
