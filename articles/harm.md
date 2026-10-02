# Testing whether a program harmed a few people

## The worry

Suppose a public benefits agency wants to make recertification stricter:
every household must document its income again, in person, within a
shorter window. The agency expects most households to clear the new
process with no trouble. Its critics worry about the few who will not —
a household that misses an appointment and loses its benefits for
months. Those households would see their food security fall sharply,
while everyone else would see no change.

The agency runs a pilot. Within each of its local offices, it assigns
half of the households at random to the new process and leaves the other
half on the old one. The question the critics care about is not whether
the new process lowered food security *on average*. It is whether the
new process harmed *anyone*. A harm confined to a few households barely
moves an average, so an analyst who compares averages can find nothing
while those households are hurt.

In this vignette we test the hypothesis that the new process changed no
household’s food security, with a test built to reject it when a few
households were harmed. The test uses scores that put their weight on
the bottom of each office’s outcomes, and it counts only evidence that
households on the stricter process ended up lower than chance would put
them.

## A simulated study

We simulate a pilot with 20 offices of 40 households each. In each
office, 20 households are assigned to the stricter process. The outcome
is a measure of food security, scaled so that the old process gives a
standard deviation of 1. Of the 400 households assigned to the stricter
process, 20 lose benefits, and their food security falls by 3. The other
380 are unaffected.

``` r

library(riposte)
set.seed(20261003)

n_office <- 20L
n_per    <- 40L
office   <- factor(rep(seq_len(n_office), each = n_per))

# assign half of each office to the stricter process, as the pilot would
assign_within_office <- function() {
  strict <- integer(n_office * n_per)
  for (o in levels(office)) {
    ix <- which(office == o)
    strict[ix][sample.int(n_per, n_per %/% 2L)] <- 1L
  }
  strict
}
strict <- assign_within_office()

# food security under the old process
food_old <- rnorm(n_office * n_per)

# 20 of the households on the stricter process lose benefits
harmed <- sample(which(strict == 1L), 20L)
food <- food_old
food[harmed] <- food[harmed] - 3

pilot <- data.frame(food = food, strict = strict, office = office)
```

Because we simulated the pilot, we know which households were harmed. An
analyst of a real pilot would not.

``` r

pilot[sort(harmed), ]
#>           food strict office
#> 34  -3.8588509      1      1
#> 35  -3.9851258      1      1
#> 103 -1.7338174      1      3
#> 149 -3.6665611      1      4
#> 168 -1.5658593      1      5
#> 176 -2.7702512      1      5
#> 230 -2.5621308      1      6
#> 282 -3.3759093      1      8
#> 347 -0.8407915      1      9
#> 379 -4.1428603      1     10
#> 422 -2.3290794      1     11
#> 509 -1.7050327      1     13
#> 525 -3.5502376      1     14
#> 543 -4.1134667      1     14
#> 583 -3.4844511      1     15
#> 617 -2.6330299      1     16
#> 690 -3.5850251      1     18
#> 699 -2.5017404      1     18
#> 744 -3.3911147      1     19
#> 750 -2.9444834      1     19
```

Average food security among households on the stricter process is -0.18,
and among households on the old process it is -0.03. The 20 harmed
households lower the average of the 400 assigned households by
$`20 \times 3 / 400 = 0.15`$, which is small next to the standard
deviation of 1.

## What a difference in means reports

We first compare averages within offices, by giving
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
the raw outcome as its only representation. A re-randomization of the
pilot is another assignment the pilot could have drawn, again putting
half of each office on the stricter process at random. The p-value is
the share of 999 re-randomizations, plus the observed assignment, whose
difference in averages is at least as far from zero as the observed one:

``` r

p_means <- riposte_test(food ~ strict | office, data = pilot,
                        representations = list(raw = function(y) y),
                        statistic = "quadratic", nresample = 999, seed = 1)$p.value
p_means
#> [1] 0.06
```

The p-value is 0.06. A researcher who sets aside the hypothesis of no
effect only at p-values of 0.05 or below would not set it aside here,
although 20 households lost 3 standard deviations of food security.

## Scores that weight the bottom of each office

Number the households in an office from the bottom: the one with the
lowest food security has rank $`k = 1`$ and the one with the highest has
rank $`k = n_b`$, the number of households in the office. The polynomial
rank score of the main vignette
([`vignette("riposte")`](https://bowers-illinois-edu.github.io/riposte/articles/riposte.md))
gives the household at rank $`k`$ the score
$`(k / (n_b + 1))^{\zeta - 1}`$, where $`\zeta`$ is a tuning value of at
least 1 that we choose. That score puts its weight on the top of the
office, and larger values of $`\zeta`$ put more of the weight on the
very top.

For harm we want the weight on the bottom.
[`riposte_poly_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_reps.md)
with `tail = "lower"` gives the household at rank $`k`$ the score

``` math
-\left(\frac{n_b + 1 - k}{n_b + 1}\right)^{\zeta - 1}.
```

Here $`n_b + 1 - k`$ is the household’s rank counted from the top.
Moving up one rank changes the score most at the bottom of the office
and hardly at all near the top. The minus sign makes the score rise with
food security, as the raw outcome, the rank, and the top-weighted score
do. So for all of these scores, lower treated outcomes mean smaller
treated-score sums, and one direction, `alternative = "less"`, describes
harm for every one of them.

We still have to choose the values of $`\zeta`$. The main vignette
derives two facts that guide the choice. First, the score at $`\zeta`$
puts most of its weight on roughly the $`n_b / (\zeta + 1)`$ households
nearest its end of the office. Second, for $`\zeta`$ = 2, 5, 14, and 41,
the values of $`2\zeta - 1`$ are 3, 9, 27, and 81, each three times the
one before. The main vignette shows that the correlation of two scores
depends only on the ratio of their values of $`2\zeta - 1`$, so values
spaced by a constant ratio give each neighbouring pair the same
correlation, and no two scores are near-copies of each other. With 40
households per office, $`\zeta`$ = 5, 14, and 41 aim at about the lowest
$`40/6 \approx 7`$, $`40/15 \approx 3`$, and $`40/42 \approx 1`$
households. We leave out $`\zeta = 2`$, the Wilcoxon score, which
responds to a shift of the whole distribution rather than to its lowest
values.

``` r

low <- riposte_poly_reps(c(5, 14, 41), tail = "lower")
names(low)
#> [1] "polylow5"  "polylow14" "polylow41"
```

## A test that counts only the harmful direction

By default,
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
counts evidence against the hypothesis of no effect in either direction.
With the bottom-weighted scores, it would find evidence when households
on the stricter process are *over*-represented at the bottom of their
offices, as harm would make them, and also when they are
*under*-represented there, as a process that lifted the worst-off
households would make them. That two-sided test gives

``` r

two_sided <- riposte_test(food ~ strict | office, data = pilot,
                          representations = low, nresample = 999, seed = 1)
two_sided$p.value
#> [1] 0.156
```

The critics’ question has one direction. With `alternative = "less"`,
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
counts only evidence that the treated-score sums are below what
re-randomization gives on average. For each value of $`\zeta`$, it sums
the scores of the households on the stricter process, after subtracting
each office’s average score. It divides each sum by its standard
deviation across re-randomizations, so that the three sums share one
scale, and takes the most negative of the three standardized sums. It
reports that most negative sum with its sign reversed as the statistic,
so a large statistic is evidence of harm. The p-value is the share of
re-randomizations, the observed assignment included, whose statistic is
at least the observed one.

``` r

one_sided <- riposte_test(food ~ strict | office, data = pilot,
                          representations = low, alternative = "less",
                          nresample = 999, seed = 1)
one_sided
#> riposte test of the sharp null of no effect
#>   800 units in 20 blocks; 999 re-randomizations
#>   combination: max
#>   alternative: less (treated-score sums below their re-randomization mean)
#>   statistic = 2.1534,  p-value = 0.0380
```

The one-sided p-value is 0.038, against 0.156 for the two-sided test and
0.06 for the difference in means. A researcher who sets aside the
hypothesis of no effect at p-values of 0.05 or below would set it aside
here and conclude that the stricter process lowered some household’s
food security. The two-sided test’s p-value is larger because it also
counts, as at least as extreme as the observed assignment,
re-randomizations in which households on the stricter process are scarce
at the bottom of their offices.

[`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
shows each score’s own one-sided p-value, which says how far up the
bottom of the offices the evidence reaches:

``` r

riposte_components(food ~ strict | office, data = pilot,
                   representations = low, alternative = "less",
                   nresample = 999, seed = 1)
#> riposte components
#>   800 units in 20 blocks; 999 re-randomizations
#>   condition number of the representation correlations: 74.8
#>   representation one-sided (less) mid-p p-values:
#>     polylow5         0.0235
#>     polylow14        0.0155
#>     polylow41        0.0365
```

## The same test when some households were helped

The p-value above tests the sharp null hypothesis that the process
changed no household’s food security. Critics might worry that the
process also helped some households, for instance those whose records
the new process corrected, and ask whether the test still means harm
when it rejects.

It does. Caughey, Dafoe, Li, and Miratrix (2023) show that one-sided
tests built from rank scores are also valid tests of the weaker
hypothesis that no unit’s outcome was lowered, which allows some
outcomes to have been raised. The reason is that the ranks in each
office are the numbers 1 through 40 whatever the process did, so the
re-randomization distribution is the same as under no effect. And when
no household was harmed, raising some outcomes on the stricter process
can only move those households up their offices’ rankings, which raises
the treated-score sums and so raises the p-value. When the one-sided
test rejects at 0.05, then, the chance of that rejection was at most
0.05 if the process harmed no one, even if it helped some. The argument
needs the default permutation engine, no covariance adjustment, and no
two households in an office with the same food security; the Details
section of
[`?riposte_test`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
gives it in full.

## How often the tests find this harm

One simulated pilot shows what one analysis looks like. To see how often
we would find the harm with each test, we repeat the pilot 100 times. We
compare the one-sided test above with the one-sided test of the raw
outcome, which asks whether the stricter process lowered average food
security. We run three versions of the pilot: the one above, one in
which 20 households are *helped* by 3 and none is harmed, and one in
which the process changes nothing.

``` r

raw <- list(raw = function(y) y)
run_pilot <- function(change) {
  strict <- assign_within_office()
  food <- rnorm(n_office * n_per)
  changed <- sample(which(strict == 1L), 20L)
  food[changed] <- food[changed] + change
  d <- data.frame(food = food, strict = strict, office = office)
  c(lower_mean = riposte_test(food ~ strict | office, data = d,
                              representations = raw, alternative = "less",
                              nresample = 499)$p.value,
    lower_tail = riposte_test(food ~ strict | office, data = d,
                              representations = low, alternative = "less",
                              nresample = 499)$p.value)
}

set.seed(3)
rejected <- sapply(c(harm_20 = -3, help_20 = 3, no_change = 0), function(change) {
  p <- replicate(100, run_pilot(change))
  rowMeans(p <= 0.05)
})
rejected
#>            harm_20 help_20 no_change
#> lower_mean    0.69    0.00      0.02
#> lower_tail    0.87    0.02      0.05
```

Each entry is the share of the 100 pilots in which the test’s p-value
was at or below 0.05. In the first column, 20 households were harmed:
with the bottom-weighted test we found the harm in 87 of 100 pilots, and
with the test of the average in 69. In the second column, 20 households
were helped and none harmed, and in the third the process changed
nothing. In both, no household was harmed, so the share estimates how
often the test reports harm when there is none, which should be at most
0.05. With 100 pilots, a share whose true value is 0.05 has a simulation
standard error of about $`\sqrt{0.05 \times 0.95 / 100} \approx
0.02`$, so shares of up to about 0.09 are consistent with a true value
of 0.05.

## Looking for a few households helped a great deal

The opposite question, whether the process helped a few households a
great deal, uses the scores that weight the top of each office,
`riposte_poly_reps(tail = "upper")`, with `alternative = "greater"`:

``` r

high <- riposte_poly_reps(c(5, 14, 41))
p_help <- riposte_test(food ~ strict | office, data = pilot,
                       representations = high, alternative = "greater",
                       nresample = 999, seed = 1)$p.value
p_help
#> [1] 0.45
```

In this pilot no household was helped, and the p-value of 0.45 gives no
reason to set aside the hypothesis of no effect.

To ask whether the process changed outcomes at either end of the
distribution, without committing to a direction, we can put the scores
for both ends into one two-sided test:

``` r

p_both <- riposte_test(food ~ strict | office, data = pilot,
                       representations = c(low, high),
                       nresample = 999, seed = 1)$p.value
p_both
#> [1] 0.268
```

## How many households, and how much

When we reject the hypothesis of no effect with the one-sided test, we
conclude that the stricter process lowered at least one household’s food
security. We learn nothing from that test about how many households were
harmed or by how much. We cannot simply count the stricter-process
households at the bottom of their offices, because some of them would
have been at the bottom under either process.

Kim, Su, Bowers, and Li, in “Randomization Tests for Distributions of
Individual Treatment Effects via Combined Rank Statistics” (arXiv
2605.08027, 2026), use the same polynomial rank scores to test
hypotheses of the form “at most $`j - 1`$ of the households on the
stricter process lost more than $`c`$.” Their `CMRSS` package computes
these tests (`remotes::install_github("bowers-illinois-edu/CMRSS")`).
Its hypotheses bound how many units had effects *above* $`c`$, so we
give it the negated outcome. A household that lost $`c`$ in food
security gained $`c`$ in the negated outcome.

``` r

n_strict <- sum(pilot$strict)
# CMRSS takes the polynomial score separately for each office
methods <- lapply(c(5, 14, 41), function(zeta)
  lapply(seq_len(n_office), function(o)
    list(name = "Polynomial", r = zeta, std = TRUE, scale = FALSE)))

# p-value for "at most j - 1 households on the stricter process lost more
# than c"; HiGHS warns about coefficients near zero that it ignores
p_at_most <- function(j, c) {
  set.seed(4)
  res <- suppressWarnings(CMRSS::pval_comb_block(
    pilot$strict, -pilot$food, k = n_strict - j + 1, c = c,
    block = pilot$office, methods.list.all = methods,
    null.max = 2000, opt.method = "ILP_highs"))
  res[["p.value"]]
}
```

With $`j = 1`$ and $`c = 0`$, the hypothesis is that no household on the
stricter process lost anything, the same hypothesis as in the section on
helped households. With offices of equal size, `CMRSS` computes the same
statistic as the one-sided
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
above, so the two p-values differ only because each package draws its
own re-randomizations. Larger $`j`$ asks whether at least $`j`$
households were harmed:

``` r

sapply(c(at_least_1 = 1, at_least_3 = 3, at_least_5 = 5), p_at_most, c = 0)
#> at_least_1 at_least_3 at_least_5 
#>     0.0330     0.0665     0.1110
```

With $`j = 1`$ and $`c`$ above zero, the hypothesis is that no household
on the stricter process lost more than $`c`$. The largest $`c`$ at which
the p-value is at or below 0.05 is a lower bound on the worst loss:

``` r

losses <- seq(0, 1, by = 0.1)
p_worst <- sapply(losses, p_at_most, j = 1)
rbind(c = losses, p = p_worst)
#>    [,1]  [,2]   [,3]   [,4]  [,5]  [,6]   [,7] [,8]  [,9] [,10]  [,11]
#> c 0.000 0.100 0.2000 0.3000 0.400 0.500 0.6000  0.7 0.800 0.900 1.0000
#> p 0.033 0.068 0.1085 0.1545 0.222 0.443 0.5215  0.8 0.885 0.951 0.9905
worst_bound <- max(c(-Inf, losses[p_worst <= 0.05]))
```

In this pilot the largest such $`c`$ on the grid is 0. So the data let a
researcher conclude, at the 0.05 level, that at least one household lost
some food security, but not that any household lost more than 0.1,
although 20 households lost 3. The tests use only ranks. A harmed
household that falls to the bottom of its office has the lowest rank
whether it lost 3 or 30, so the ranks carry little information about the
size of a loss, and a test built on them can show that someone was
harmed long before it can show by how much.

## References

Caughey, D., Dafoe, A., Li, X., and Miratrix, L. (2023). Randomisation
inference beyond the sharp null: bounded null hypotheses and quantiles
of individual treatment effects. *Journal of the Royal Statistical
Society Series B: Statistical Methodology*, 85(5), 1471–1491.
<https://doi.org/10.1093/jrsssb/qkad080>

Kim, D., Su, Y., Bowers, J., and Li, X. (2026). Randomization tests for
distributions of individual treatment effects via combined rank
statistics. arXiv 2605.08027.
