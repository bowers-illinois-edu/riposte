# Why the large-sample Cauchy combination is truncated at 0.9

## A combined p-value of 1 next to a difference in means of 0.028

Wallsten and Nteta (2016) ran a survey experiment in which 36
respondents, 18 of them treated, answered on a 5-point scale coded 0,
0.25, 0.5, 0.75, and 1. Sarkar and Coppock (2026) reanalyzed it with
other experiments on religious messages, and the cleaned data in their
replication archive give these counts:

``` r

vals <- c(0, 0.25, 0.5, 0.75, 1)
wal <- data.frame(Y = c(rep(vals, c(10, 2, 3, 2, 1)), rep(vals, c(5, 2, 0, 8, 3))),
                  Z = rep(0:1, each = 18))
table(answer = wal$Y, treated = wal$Z)
#>       treated
#> answer  0  1
#>   0    10  5
#>   0.25  2  2
#>   0.5   3  0
#>   0.75  2  8
#>   1     1  3
```

The difference in means gives a p-value of 0.028. riposte’s Cauchy
combination, with its p-value from the large-sample formula that riposte
used before version 0.0.0.9008, gives exactly 1:

``` r

t.test(Y ~ Z, data = wal)$p.value
#> [1] 0.02825938
riposte_test(Y ~ Z, data = wal, statistic = "cauchy", engine = "asymptotic",
             cauchy_truncation = 1)$p.value
#> [1] 1
```

A combined p-value of 1 says that the observed result is the least
remarkable result any assignment of treatment could have produced. This
vignette shows where that 1 comes from, how the truncated conversion of
Gui, Jiang, and Wang (2025) removes it, and why riposte uses their
truncation level of 0.9 by default. We will show that one of the six
representations has a correct p-value of exactly 1, that the untruncated
Cauchy combination turns any p-value of 1 into a combined p-value of 1,
and that a smaller truncation level than 0.9 would trade a little more
power for a few more false positives.

## The untruncated Cauchy combination lets one p-value of 1 decide

The Cauchy combination of Liu and Xie (2020) converts each
representation’s p-value $`p`$ to $`\tan((0.5 - p)\pi)`$, averages the
converted values, and turns the average back into a p-value with the
standard Cauchy distribution, for which a variable $`C`$ exceeds $`x`$
with probability $`P(C > x) = 1/2 - \arctan(x)/\pi`$. The conversion
makes a small p-value into a large number, so that one representation
with strong evidence is not diluted by others with none:

``` r

p <- c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999, 1)
data.frame(p = p, converted = format(tan((0.5 - p) * pi), digits = 4))
#>       p  converted
#> 1 0.001  3.183e+02
#> 2 0.010  3.182e+01
#> 3 0.050  6.314e+00
#> 4 0.500  0.000e+00
#> 5 0.950 -6.314e+00
#> 6 0.990 -3.182e+01
#> 7 0.999 -3.183e+02
#> 8 1.000 -1.633e+16
```

The conversion treats p-values near 1 the same way with the sign
reversed: 0.999 becomes about $`-318`$, and 1 becomes minus infinity (R
prints a huge finite number because its value of $`\pi/2`$ is rounded).
One converted value of minus infinity makes the average minus infinity,
and the combined p-value is then 1 whatever the other p-values are.
riposte’s default representations are six: the raw outcome, its rank,
mean distance, mean rank distance, max distance, and Huber (the outcome
minus its median, divided by the median of the absolute differences from
the median, and capped at plus or minus 1.345). These are their six
large-sample p-values in the Wallsten and Nteta data:

``` r

riposte_test(Y ~ Z, data = wal, statistic = "cauchy", engine = "asymptotic")$component_p
#>            raw           rank      mean_dist mean_rank_dist       max_dist 
#>     0.03030884     0.02925259     0.11732820     0.18478113     1.00000000 
#>          huber 
#>     0.03011563
```

The max distance’s p-value is exactly 1. The max distance scores each
answer by its distance to the farther end of the observed range, which
measures how far the answer is from the middle of the scale. To compare
the treated respondents’ average score with everyone’s, we need each
answer’s score:

| answer       | 0   | 0.25 | 0.5 | 0.75 | 1   |
|--------------|-----|------|-----|------|-----|
| max distance | 1   | 0.75 | 0.5 | 0.75 | 1   |

The 18 treated respondents include 8 at an end and 10 one step in, so
their average score is $`(8 \times 1 + 10 \times 0.75)/18 = 31/36`$. All
36 include 19 at an end, 14 one step in, and 3 in the middle, so their
average is $`(19 \times 1 + 14 \times 0.75 + 3 \times 0.5)/36 = 31/36`$.
The treated average equals the overall average, so the standardized
difference is exactly 0, and its two-sided p-value, the chance of a
difference at least as far from 0, is 1. That p-value is correct: on the
max distance the arms do not differ at all. On outcomes with few values
and few units, a treated average equals the overall average exactly in a
noticeable share of assignments, and this can happen to any
representation.

## The truncated conversion of Gui, Jiang, and Wang

Gui, Jiang, and Wang (2025, Section 2.2) convert each p-value with the
quantile function of a Cauchy distribution cut off below. They keep a
share $`t`$ of the Cauchy distribution, above the point
$`c = \tan((0.5 - t)\pi)`$. Write $`F_C`$ for the Cauchy distribution
function, so $`F_C(c) = 1 - t`$. The cut-off distribution has
distribution function $`F(x) = (F_C(x) - (1 - t))/t`$ for $`x \ge c`$.
Its quantile for $`p`$ is the $`x`$ with $`F(x) = 1 - p`$, that is
$`F_C(x) = 1 - tp`$, and since the Cauchy quantile of $`1 - q`$ is
$`\tan((0.5 - q)\pi)`$, the conversion is $`x = \tan((0.5 - tp)\pi)`$.
With $`t = 0.9`$, a p-value of 1 converts to $`\tan(-0.4\pi) = -3.08`$
instead of minus infinity. The cut-off distribution’s upper tail is the
Cauchy’s divided by $`t`$. With $`n`$ p-values whose converted values
sum to $`S`$, Gui, Jiang, and Wang take the combined p-value to be
$`\min(1, n P(C > S)/t)`$, with $`C`$ a standard Cauchy variable (their
Section 2.3). The reasoning is that for variables with heavy upper
tails, a large sum almost always comes from one large term, so the
chance that the sum of $`n`$ converted values exceeds $`S`$ is about
$`n`$ times the chance that one of them does, and one converted value
exceeds $`S`$ with probability $`P(C > S)/t`$. With one p-value this
returns that p-value exactly.
[`riposte_truncated_cauchy()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_truncated_cauchy.md)
computes it, and
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
uses it for the large-sample Cauchy combination and hybrid:

``` r

riposte_test(Y ~ Z, data = wal, statistic = "cauchy", engine = "asymptotic")$p.value
#> [1] 0.05736861
```

The max distance’s p-value of 1 now converts to $`-3.08`$ instead of
minus infinity, and the raw outcome, rank, and Huber p-values near 0.03
still count. Fang, Chang, Park, and Tseng (2023) proposed an earlier
truncation of the same kind, replacing every p-value above 0.99 with
0.99, and noted that p-values near 1 are common with discrete data.

## What the truncation level changes

The choice of $`t`$ is easiest to see in the rule for rejecting. With
six p-values, the combined p-value is at most 0.05 exactly when
$`P(C > S) \le 0.05t/6`$, that is, when
$`S \ge \cot(\pi \cdot 0.05t/6)`$, since $`P(C > x) = q`$ exactly when
$`x = \cot(\pi q)`$. Since each converted value is $`\cot(\pi t p_k)`$,
multiplying both sides by $`\pi t`$ gives
$`\sum_k \pi t \cot(\pi t p_k) \ge \pi t \cot(\pi t \cdot 0.05/6)`$. We
call $`\pi t \cot(\pi t p)`$ the contribution of a p-value $`p`$. The
rule is then: reject when the six contributions add up to at least the
contribution of a single p-value of $`0.05/6`$. Because
$`\pi t \cot(\pi t p)`$ is close to $`1/p`$ for small $`p`$, that
threshold is close to $`6/0.05 = 120`$ for any $`t`$.

``` r

contribution <- function(p, t) pi * t / tan(pi * t * p)
p <- c(0.01, 0.05, 0.5, 0.8, 1)
round(cbind(p, `t = 0.9` = contribution(p, 0.9), `t = 0.8` = contribution(p, 0.8),
            `t = 0.5` = contribution(p, 0.5), `1/p` = 1 / p), 2)
#>         p t = 0.9 t = 0.8 t = 0.5    1/p
#> [1,] 0.01   99.97   99.98   99.99 100.00
#> [2,] 0.05   19.87   19.89   19.96  20.00
#> [3,] 0.50    0.45    0.82    1.57   2.00
#> [4,] 0.80   -2.34   -1.18    0.51   1.25
#> [5,] 1.00   -8.70   -3.46    0.00   1.00
```

A small p-value contributes about $`1/p`$ whatever $`t`$ is, so $`t`$
changes only what large p-values contribute. At $`t = 0.9`$ a p-value of
1 subtracts 8.7, at $`t = 0.8`$ it subtracts 3.5, and at $`t = 0.5`$ it
subtracts nothing. As $`t`$ approaches 0 every p-value contributes
$`1/p`$, and the rule becomes $`\sum_k 1/p_k \ge 6/0.05`$, that is,
$`6/\sum_k (1/p_k) \le 0.05`$. That is the harmonic mean p-value of
Wilson (2019): the number of p-values divided by the sum of their
reciprocals. Lowering $`t`$ only raises the contributions of large
p-values, so a smaller $`t`$ rejects in every experiment that a larger
$`t`$ rejects, and in some more. Power and the false positive rate
therefore rise together as $`t`$ falls, and the smallest p-value over
several values of $`t`$ is always the p-value at the smallest of them.

## How much power, and how many false positives

`dev/truncation-sims.R` in the package source measures both. For false
positives it holds outcomes fixed, so that the treatment changes
nothing, and reassigns treatment 20,000 times: once for the Wallsten and
Nteta data and once for 40 outcomes drawn from a skewed (lognormal)
distribution, half treated. For power it draws 4,000 experiments with 40
lognormal outcomes, half treated, in which the treatment either adds
0.55 to every treated outcome or moves treated outcomes 2.1 times as far
from the control median. Every rule is applied to the same six
large-sample p-values. Each entry is the share of experiments rejected
at 0.05:

|  | false positives, Wallsten | false positives, lognormal | power, shift | power, spread |
|:---|---:|---:|---:|---:|
| untruncated (Liu and Xie) | 0.043 | 0.042 | 0.474 | 0.542 |
| t = 0.95 | 0.047 | 0.043 | 0.484 | 0.554 |
| t = 0.9 | 0.047 | 0.044 | 0.486 | 0.558 |
| t = 0.8 | 0.049 | 0.044 | 0.490 | 0.561 |
| t = 0.7 | 0.050 | 0.045 | 0.493 | 0.564 |
| t = 0.5 | 0.051 | 0.045 | 0.496 | 0.566 |
| harmonic mean (t near 0) | 0.053 | 0.046 | 0.499 | 0.569 |

The standard error of a rate near 0.05 from 20,000 reassignments is
$`\sqrt{0.05 \times 0.95 / 20000} = 0.0015`$, and of a power near 0.5
from 4,000 experiments about 0.008. The six p-values in these designs
are correlated. With six independent p-values, the case in which these
rules exceed 0.05 by the most in Gui, Jiang, and Wang’s Figure 2, one
million simulated sets give these false positive rates:

| rule                      | false positive rate |
|:--------------------------|--------------------:|
| untruncated (Liu and Xie) |              0.0498 |
| t = 0.95                  |              0.0536 |
| t = 0.9                   |              0.0551 |
| t = 0.8                   |              0.0569 |
| t = 0.7                   |              0.0581 |
| t = 0.5                   |              0.0595 |
| harmonic mean (t near 0)  |              0.0608 |

The untruncated combination rejects too rarely in the Wallsten data,
because a p-value of exactly 1 occurs in many reassignments and forces
the combined p-value to 1 each time. Truncating at 0.9 brings the rate
close to 0.05 and raises power. Moving from 0.9 to smaller levels adds
about a point of power at most, while the false positive rate keeps
rising. With independent p-values it is already 0.055 at 0.9 and 0.059
at 0.5. Gui, Jiang, and Wang recommend 0.9 from their own simulations
(their Section 6), and riposte uses it.

## Several truncation levels at once

Since the smallest p-value over several truncation levels is the p-value
at the smallest level, choosing the best level after seeing the data is
the same as always using the smallest, which has the most false
positives above. Taking the p-value from re-randomization removes that
problem. When the treatment changed no one’s outcome, the observed
assignment is one more random draw from the assignment process, so its
statistic is no more likely than any reassignment’s to be the most
extreme. A p-value that counts the share of reassignments at least as
extreme as the observed statistic is then at most 0.05 no more than 5%
of the time, whatever the statistic is. So the simulation also took each
statistic’s p-value from 999 reassignments, in 2,000 experiments per
column:

|  | false positives, Wallsten | power, shift | power, spread |
|:---|---:|---:|---:|
| t = 0.9, large-sample | 0.050 | 0.488 | 0.584 |
| smallest over t = 0.5, 0.7, 0.9, large-sample | 0.056 | 0.498 | 0.591 |
| t = 0.9, re-randomization | 0.056 | 0.516 | 0.612 |
| smallest over t = 0.5, 0.7, 0.9, re-randomization | 0.056 | 0.515 | 0.613 |

The standard errors here are about 0.005 for the false positive rate and
0.011 for power. Under re-randomization, one level and several levels
reject at nearly the same rates, so using several levels gains nothing.
Taking the p-value from re-randomization rather than from the
large-sample formula raised power by about 3 percentage points at
$`t = 0.9`$.

## The re-randomization Cauchy combination

With the default `engine = "permute"`, riposte gives each representation
a mid-p value: the share of reassignments whose statistic is more
extreme than the observed one, counting reassignments tied with it as
one half. A mid-p value is never exactly 1. riposte averages the
converted mid-p values and takes the share of reassignments whose
average is at least the observed one. For the reason given in the
previous section, the test’s false positive rate is then at most the
nominal level whatever the conversion. `cauchy_truncation` therefore
applies only to `engine = "asymptotic"`.

## References

Fang, Y., Chang, C., Park, Y., and Tseng, G. C. (2023). Heavy-tailed
distribution for combining dependent p-values with asymptotic
robustness. *Statistica Sinica*, 33, 1115-1142.
<doi:10.5705/ss.202022.0046>

Gui, L., Jiang, Y., and Wang, J. (2025). Aggregating dependent signals
with heavy-tailed combination tests. *Biometrika*, 112(4), asaf038.
<doi:10.1093/biomet/asaf038>

Liu, Y., and Xie, J. (2020). Cauchy combination test: a powerful test
with analytic p-value calculation under arbitrary dependency structures.
*Journal of the American Statistical Association*, 115(529), 393-402.
<doi:10.1080/01621459.2018.1554485>

Long, M., Li, Z., Zhang, W., and Li, Q. (2023). The Cauchy combination
test under arbitrary dependence structures. *The American Statistician*,
77(2), 134-142. <doi:10.1080/00031305.2022.2116109>

Sarkar, R., and Coppock, A. (2026). The effects of religious messages
and endorsements on political attitudes: a meta-reanalysis. *American
Political Science Review*, advance online publication.
<doi:10.1017/S0003055426101695>. Replication archive:
<doi:10.7910/DVN/NHDFFS>

Wallsten, K., and Nteta, T. M. (2016). For you were strangers in the
land of Egypt: clergy, religiosity, and public opinion toward
immigration reform in the United States. *Politics and Religion*, 9(3),
566-604. <doi:10.1017/S1755048316000444>

Wilson, D. J. (2019). The harmonic mean p-value for combining dependent
tests. *Proceedings of the National Academy of Sciences*, 116(4),
1195-1200. <doi:10.1073/pnas.1814092116>
