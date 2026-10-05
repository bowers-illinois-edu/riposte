# Detecting effects a difference in means can miss

## The problem

Suppose you run a randomized program evaluation. You randomize within
blocks — schools, clinics, counties — so that comparable units are
compared, and you want to know whether the program did anything. The
reflex is to compare the average outcome in the treated group to the
average in the control group. That comparison answers one question: did
the program shift the *average*?

A program can change outcomes in ways the average hides. It might help a
few people a great deal and leave everyone else alone. It might raise
some places and lower others, so the average barely moves. It might
widen the spread — making good outcomes better and bad outcomes worse —
without changing the center. In each case a difference in means reports
“no effect,” and you would conclude the program did nothing when it did
something.

`riposte` tests for an effect by looking at the outcome several ways at
once — its values, its ranks, how far apart units are, a bounded version
— and combining what each view sees. By default it reports a single
p-value that is *exact*: it comes from the actual randomization you ran,
not from a large-sample approximation. So if the program changed no
one’s outcome, the chance that the p-value falls at or below 0.05 is at
most 0.05, however many units or blocks you have.

## A worked example

Here is a block-randomized experiment where the program changes the
*spread* of the outcome but not its average. Treated units have the same
mean as controls and a larger variance.

``` r

library(riposte)

B  <- 12L   # blocks
nb <- 30L   # units per block
block <- factor(rep(seq_len(B), each = nb))

# assign half of each block to treatment, at random, the way the study would
z <- integer(B * nb)
for (b in levels(block)) {
  ix <- which(block == b)
  z[ix][sample.int(nb, nb %/% 2L)] <- 1L
}

# treatment changes the spread, not the mean
y0 <- rnorm(B * nb, mean = 0, sd = 1)
y1 <- rnorm(B * nb, mean = 0, sd = 1.8)
y  <- ifelse(z == 1, y1, y0)

dat <- data.frame(outcome = y, treated = z, block = block)
```

A difference in means sees almost nothing here, because the averages are
the same. `riposte` looks at the whole distribution:

``` r

riposte_test(outcome ~ treated | block, data = dat,
             statistic = "screen", nresample = 999, seed = 1)
#> riposte test of the sharp null of no effect
#>   360 units in 12 blocks; 999 re-randomizations
#>   screen (shrink, lambda = 0.017); condition number 235.0
#>   combination: quadratic
#>   statistic = 41.5062,  p-value = 0.0010
```

The test reports a small p-value: the program *did* change the outcome,
through its spread.

## Choosing how to calculate the p-value

The formula specifies the assignment design. `outcome ~ treated | block`
means that treatment was randomized within blocks, with each block’s
treated count fixed. `outcome ~ treated` means complete randomization,
with only the total treated count fixed. Either design can be used with
either calculation below.

With `engine = "permute"`, the default, we generate assignments under
that design and compare the observed statistic with the statistics at
those assignments. With `engine = "asymptotic"`, we instead use an
approximation: a chi-square distribution for the quadratic statistic, or
a standard Cauchy distribution for the Cauchy statistic.

``` r

riposte_test(outcome ~ treated | block, data = dat,
             statistic = "quadratic", engine = "asymptotic")
#> riposte test of the sharp null of no effect
#>   360 units in 12 blocks
#>   asymptotic engine (chi-square, df = 6; no re-randomization)
#>   combination: quadratic
#>   statistic = 41.5313,  approximate p-value = 2.27523e-07
riposte_test(outcome ~ treated | block, data = dat,
             statistic = "cauchy", engine = "asymptotic")
#> riposte test of the sharp null of no effect
#>   360 units in 12 blocks
#>   asymptotic engine (standard Cauchy; no re-randomization)
#>   combination: cauchy
#>   statistic = 372593114.2778,  approximate p-value = 9.49233e-10
```

For the quadratic approximation, the degrees of freedom are the number
of independent score sums retained in the covariance inverse. Repeating
a score does not add another degree of freedom. For the Cauchy
calculation, we first obtain each representation’s p-value from its
squared standardized statistic and a chi-square distribution on one
degree of freedom. We then average the Cauchy transforms of those
p-values and calculate the standard Cauchy upper tail. The two combined
p-values need not equal their permutation counterparts.

The asymptotic engine uses no random draws, so `seed` and `nresample`
have no effect. Its printed result identifies the approximation. It
supports the unadjusted quadratic and Cauchy tests; the max test, the
default screen, and covariance adjustment require other methods. The
approximations do not give the permutation test’s finite-sample level
guarantee. In particular, individual scores that dominate the variance
can make the normal approximation behind the chi-square calculation
inaccurate. Dependence among the individual p-values also means that the
Cauchy reference is an approximation.

### The Hybrid Cauchy example

The paper’s Hybrid test combines seven p-values: one for each of the six
representations, and one from their quadratic combination. We give each
p-value weight `1 / 7` by averaging its Cauchy transform with the other
six. The existing functions supply all seven inputs:

``` r

six <- riposte_test(outcome ~ treated | block, data = dat,
                    statistic = "cauchy", engine = "asymptotic")
quad <- riposte_test(outcome ~ treated | block, data = dat,
                     statistic = "quadratic", engine = "asymptotic")
hybrid_inputs <- c(six$component_p, quadratic = quad$p.value)
hybrid_statistic <- mean(riposte_acat_term(hybrid_inputs))
hybrid_p <- pcauchy(hybrid_statistic, lower.tail = FALSE)
hybrid_p
#> [1] 1.106669e-09
```

`six$component_p` contains the six individual p-values. Using
`six$p.value` instead would combine two p-values at equal weights and
give the quadratic test half the weight. The example uses the standard
Cauchy tail to approximate the combined p-value; it does not compare the
Hybrid statistic with its permutation distribution. As with the
six-component approximation, an input of exactly one gives a negative
infinite term. Opposing infinite terms leave the sum undefined.

## What each view of the outcome contributes

[`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
shows the evidence each representation carries on its own, so you can
see *which way* the effect shows up. These diagnostics use permutation
mid-p values, not the asymptotic p-values used in the example above:

``` r

riposte_components(outcome ~ treated | block, data = dat,
                   nresample = 999, seed = 1)
#> riposte components
#>   360 units in 12 blocks; 999 re-randomizations
#>   condition number of the representation correlations: 235.0
#>   representation two-sided mid-p p-values:
#>     raw              0.6835
#>     rank             0.9375
#>     mean_dist        0.0005
#>     mean_rank_dist   0.0005
#>     max_dist         0.0005
#>     huber            0.9135
```

The distance representations (how far each unit sits from the others)
light up, because a spread change moves units apart; the raw values and
ranks do not, because the average and the ordering are unchanged. A
difference in means uses only the first of these.

## Combining the views: the quadratic and the Cauchy

There are two natural ways to combine several views, and they suit
different situations.

- The **quadratic** combination uses how the views move together. When
  an effect leaves a weak trace in several correlated views, the
  quadratic adds those traces up and finds it. This is the default’s
  workhorse.
- The **Cauchy** combination uses each view only on its own. It wins
  when one view carries a strong, clear signal and the others are noise.

You rarely have to choose by hand. The **screen** (the default) reads
how much the views overlap — a quantity fixed by your data and design,
not by the coin flips of the randomization — and leans toward the
quadratic when the views are well behaved and toward the Cauchy when
they overlap so much that combining them would be unstable. Because that
quantity does not depend on the realized assignment, choosing from it
does not disturb the test’s exactness.

``` r

# force a single combination if you want one
riposte_test(outcome ~ treated | block, data = dat,
             statistic = "cauchy", nresample = 999, seed = 1)$p.value
#> [1] 0.001
```

## Rank scores that weight the top of each block

Some programs help a few people a great deal and leave the rest alone. A
rank score that gives the top of each block’s ranking more weight than
the rest responds to that pattern through the ranks alone. `riposte`
offers the polynomial rank score of Kim, Su, Bowers, and Li (arXiv
2605.08027): the unit at rank $`k`$ in a block of $`n_b`$ units gets

``` math
\left(\frac{k}{n_b + 1}\right)^{\zeta - 1}.
```

At $`\zeta = 2`$ the score is the rank divided by $`n_b + 1`$, which
gives the Wilcoxon rank-sum test. As $`\zeta`$ grows, more of each
block’s total score sits on its highest-ranked units. The score always
lies between 0 and 1, so blocks of different sizes enter on the same
scale.

[`riposte_poly_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_reps.md)
returns one representation per value of $`\zeta`$. Its default,
$`\zeta`$ = 2, 7, 12, 17, 22, is the set in Bowers and Burton’s
rank-score tables. Here are the five scores on the spread-change
experiment above, combined and one at a time:

``` r

riposte_test(outcome ~ treated | block, data = dat,
             representations = riposte_poly_reps(), nresample = 999, seed = 1)
#> riposte test of the sharp null of no effect
#>   360 units in 12 blocks; 999 re-randomizations
#>   screen (shrink, lambda = 0.042); condition number 134869.4
#>   combination: quadratic
#>   statistic = 21.9384,  p-value = 0.0010
riposte_components(outcome ~ treated | block, data = dat,
                   representations = riposte_poly_reps(), nresample = 999, seed = 1)
#> riposte components
#>   360 units in 12 blocks; 999 re-randomizations
#>   condition number of the representation correlations: 134869.4
#>   representation two-sided mid-p p-values:
#>     poly2            0.9375
#>     poly7            0.0045
#>     poly12           0.0025
#>     poly17           0.0015
#>     poly22           0.0025
```

The score at $`\zeta = 2`$ detects a shift in location, and the two
groups have the same mean, so it finds nothing. The program widens the
spread, which puts more treated units at the top of each block, so the
scores at larger $`\zeta`$ find it.

### Choosing your own values of $`\zeta`$

Pass any numbers of at least 1 to
[`riposte_poly_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_reps.md).
They need not be whole numbers. Three questions decide which numbers to
pass.

The smallest value is 2. At $`\zeta = 1`$ every unit scores 1, and
`riposte` drops a score that is constant within every block.

The largest value depends on how few units per block you want the test
to be able to find gains in. Write $`x = k/(n_b + 1)`$ for a unit’s rank
as a fraction of its block. If $`x`$ were spread evenly between 0 and 1,
the score-weighted average distance from the top of the block would be
$`1/(\zeta + 1)`$ of the block, so the score at $`\zeta`$ concentrates
on roughly the top $`n_b/(\zeta + 1)`$ units. To find gains confined to
about $`m`$ units per block, make the largest $`\zeta`$ about
$`n_b/m - 1`$. The code below computes the exact average for blocks of
30, the size in our example:

``` r

nb <- 30
k <- seq_len(nb)
zetas <- c(2, 5, 7, 12, 14, 22, 29)
where <- sapply(zetas, function(zeta) {
  s <- riposte_poly_scores(k, zeta)
  sum((nb + 1 - k) * s) / sum(s)      # score-weighted ranks from the top
})
round(rbind(zeta = zetas, exact = where, approx = nb / (zetas + 1)), 1)
#>        [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> zeta    2.0  5.0  7.0 12.0 14.0 22.0 29.0
#> exact  10.7  5.6  4.3  2.9  2.6  1.9  1.6
#> approx 10.0  5.0  3.8  2.3  2.0  1.3  1.0
```

Here `exact` counts the top unit as 1 rank from the top, and it runs 0.5
to 0.7 units above `approx`.

The values in between should not be near-copies of one another. Two
scores that move together across re-randomizations carry almost the same
information, yet the quadratic combination counts each as a separate
direction and has to invert their nearly singular correlation matrix.
With $`x`$ spread evenly, the scores at $`\zeta`$ and $`\zeta'`$
correlate at

``` math
\frac{2\sqrt{uv}}{u + v}, \qquad u = 2\zeta - 1,\; v = 2\zeta' - 1,
```

which depends only on the ratio $`v/u`$. Equal steps in $`\zeta`$ give
neighbours whose correlation climbs toward 1. Equal ratios in
$`2\zeta - 1`$ give every neighbouring pair the same correlation:
tripling $`2\zeta - 1`$ from 3 gives $`\zeta`$ = 2, 5, 14, 41, and each
neighbouring pair correlates at $`2\sqrt{3}/4 = 0.87`$. For blocks of
30, the set 2, 5, 14 reaches about 2.6 units from the top, by the
`exact` row of the table above. Within one block, the correlation of two
treated-unit sums across re-randomizations equals the correlation of the
two score vectors, so you can check any candidate set for your own block
size before running a test:

``` r

score_cor <- function(zetas, nb) {
  k <- seq_len(nb)
  S <- sapply(zetas, function(zeta) riposte_poly_scores(k, zeta))
  dimnames(S) <- list(NULL, paste0("zeta", zetas))
  cor(S)
}
round(score_cor(c(2, 7, 12, 17, 22), nb = 30), 3)
#>        zeta2 zeta7 zeta12 zeta17 zeta22
#> zeta2  1.000 0.786  0.646  0.562  0.506
#> zeta7  0.786 1.000  0.962  0.905  0.853
#> zeta12 0.646 0.962  1.000  0.985  0.957
#> zeta17 0.562 0.905  0.985  1.000  0.992
#> zeta22 0.506 0.853  0.957  0.992  1.000
round(score_cor(c(2, 5, 14), nb = 30), 3)
#>        zeta2 zeta5 zeta14
#> zeta2  1.000 0.870  0.608
#> zeta5  0.870 1.000  0.869
#> zeta14 0.608 0.869  1.000

# ratio of the largest to the smallest eigenvalue: how unstable the inverse is
condition <- function(R) { ev <- eigen(R)$values; max(ev) / min(ev) }
condition(score_cor(c(2, 7, 12, 17, 22), nb = 30))
#> [1] 134869.4
condition(score_cor(c(2, 5, 14), nb = 30))
#> [1] 69.84958
```

From $`\zeta = 7`$ up, the default set’s neighbours correlate at 0.96 or
more, and the ratio of its largest to its smallest eigenvalue is in the
hundreds of thousands. The default screen still runs on it, because it
shrinks the correlation matrix toward the identity until that ratio is
100 (the `threshold` of
[`riposte_screen_control()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_screen_control.md)).
The set 2, 5, 14 needs no shrinking. A set chosen this way, for these
blocks of 30:

``` r

riposte_test(outcome ~ treated | block, data = dat,
             representations = riposte_poly_reps(c(2, 5, 14)),
             nresample = 999, seed = 1)
#> riposte test of the sharp null of no effect
#>   360 units in 12 blocks; 999 re-randomizations
#>   screen (shrink, lambda = 0.000); condition number 69.8
#>   combination: quadratic
#>   statistic = 28.5797,  p-value = 0.0010
```

These correlations describe the scores under the null hypothesis; they
do not say which set is more powerful against a particular effect. A
larger $`\zeta`$ can add power against gains among very few units even
when it correlates highly with its neighbour under the null.

When blocks differ in size, one value of $`\zeta`$ concentrates on the
same fraction of every block, about $`1/(\zeta + 1)`$, not on the same
number of units.

### Adding rank scores to the six defaults

The six default representations already include the rank, and once each
score is centred within its block, the score at $`\zeta = 2`$ is the
rank divided by $`n_b + 1`$. When every block has the same size, the two
are exactly proportional. So leave $`\zeta = 2`$ out when you add
polynomial scores to the six:

``` r

riposte_test(outcome ~ treated | block, data = dat,
             representations = c(riposte_reps_default(),
                                 riposte_poly_reps(c(5, 14))),
             nresample = 999, seed = 1)
#> riposte test of the sharp null of no effect
#>   360 units in 12 blocks; 999 re-randomizations
#>   screen (shrink, lambda = 0.042); condition number 2386.5
#>   combination: quadratic
#>   statistic = 41.8322,  p-value = 0.0010
```

To look for a few people *harmed* by a program, use
`riposte_poly_reps(tail = "lower")`, which weights the bottom of each
block, with `riposte_test(alternative = "less")`, which counts only
evidence that treated outcomes are lower.
[`vignette("harm")`](https://bowers-illinois-edu.github.io/riposte/articles/harm.md)
works through that case.

These tests, like every test in `riposte`, test the sharp null
hypothesis that the program changed no one’s outcome. Kim, Su, Bowers,
and Li use the same scores to test a different hypothesis, about how
many units had an effect larger than a given amount. `riposte` does not
compute those tests; the `CMRSS` package does.

## Bringing in covariates

If you measured covariates that predict the outcome, you can use them to
sharpen the test. `riposte` fits a model on the **control units only**,
subtracts its prediction, and tests what is left. Fitting on controls
only means the model can never absorb the treatment effect, and
refitting it inside every re-randomization keeps the test exact for any
model you choose — ridge regression by default, but a random forest or
anything else works the same way.

``` r

dat$x1 <- rnorm(nrow(dat)); dat$x2 <- rnorm(nrow(dat))
riposte_test(outcome ~ treated | block, data = dat,
             adjust = ~ x1 + x2, nresample = 999, seed = 1)$p.value
#> [1] 0.001
```

When the covariates predict the outcome, the model removes variation in
the outcome that the treatment did not cause, which makes the test more
powerful. When they do not, the model has little to remove. In 100
simulated experiments shaped like the one above, but with treatment
raising the standard deviation from 1 to 1.35 and with two covariates of
pure noise, the unadjusted test rejected at the 0.05 level in 85 percent
of experiments and the adjusted test in 83 percent. That simulation is
one of the package’s tests, in `tests/testthat/test-power.R` in the
source repository.

Whether or not the covariates predict the outcome, the model is refit on
each re-randomization, so the chance of a p-value at or below 0.05 when
the program changed no one’s outcome stays at most 0.05.

## Cluster-randomized designs

If you assigned treatment to whole groups — all classrooms in a school,
all clinics in a county — then the groups, not the individuals, are what
was randomized. Tell `riposte` which units share a cluster and it
permutes whole clusters, so the p-value reflects the design you actually
ran. Treating the individuals as if they were independently randomized
would overstate your certainty.

``` r

riposte_test(outcome ~ treated | block, data = dat,
             clusters = "cluster_id", nresample = 999, seed = 1)
```

## What to remember

A difference in means asks one narrow question. `riposte` asks a broader
one — did the program change the distribution of outcomes, in any of
several ways — and by default answers it with a p-value that is exact
under your randomization. Reach for it when an effect might be
concentrated, offsetting, or in the spread rather than the average, and
when you want the validity of a randomization test rather than a
large-sample approximation.
