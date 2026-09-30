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
not from a large-sample approximation, so it is valid however many units
or blocks you have.

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
#>   statistic = 335333801.9563,  approximate p-value = 9.49233e-10
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

When the covariates predict the outcome, adjustment makes the test more
powerful; when they do not, it costs almost nothing, and either way the
p-value stays valid because the model is refit on each re-randomization.

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
