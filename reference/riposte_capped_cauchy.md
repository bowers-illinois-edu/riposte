# Capped Cauchy combination of p-values (the paper's rule)

Replaces every p-value above `cap` by `cap` and then combines the
p-values by Liu and Xie's Cauchy combination: the mean of \\\tan((0.5 -
p)\pi)\\, referred to the standard Cauchy distribution. This is the rule
of Bowers and Burton's code (footnote in Section 4.5 of "A More Powerful
Test for Randomized Experiments"), with `cap = 0.999`. Without a cap one
p-value of 1 converts to minus infinity and makes the combined p-value 1
whatever the others are. With the cap it converts to \\\tan(-0.499\pi) =
-318\\. riposte's default for the large-sample Cauchy and hybrid is
instead the truncation of Gui, Jiang, and Wang, which converts a p-value
of 1 to -3.08
([`riposte_truncated_cauchy()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_truncated_cauchy.md));
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
applies this cap when given `cauchy_cap`.

## Usage

``` r
riposte_capped_cauchy(p, cap = 0.999)
```

## Arguments

- p:

  numeric vector of p-values in \[0, 1\].

- cap:

  the largest p-value that enters the combination, in (0.5, 1\]; 0.999
  by default. `cap = 1` gives Liu and Xie's combination unchanged.

## Value

the combined p-value.

## Examples

``` r
# the seven p-values of the paper's example: uncapped the combination is 1
p <- c(0.0003, 0.15, 0.57, 0.05, 0.81, 0.25, 1)
riposte_capped_cauchy(p)            # 0.0030
#> [1] 0.002969598
riposte_capped_cauchy(p, cap = 1)   # 1
#> [1] 1
riposte_truncated_cauchy(p)         # riposte's default rule: 0.0021
#> [1] 0.002088631
```
