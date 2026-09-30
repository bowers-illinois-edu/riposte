# Huber's psi of a block's outcomes, standardised by the block's MAD

A bounded representation: each outcome's departure from the block
median, in units of the block's median absolute deviation (MAD, scaled
by 1.4826 as [`stats::mad()`](https://rdrr.io/r/stats/mad.html) does, so
that it equals the standard deviation for normal outcomes), held within
`[-k, k]`. With `k = 1.345`, Huber's M-estimate of location is 95
percent as efficient as the mean when the outcomes are normal.

## Usage

``` r
riposte_huber(y, k = 1.345)
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

- k:

  the cap, in MAD units. Huber's 1.345 by default.

## Value

numeric vector the same length as `y`.

## Details

Standardising within block keeps the representation informative whatever
the outcome's location and units, and the cap keeps one extreme unit
from dominating the block. When more than half of a block's outcomes
tie, the MAD is zero; the departures are then divided by the mean
absolute departure from the median instead. A block with no spread at
all gets zeros.

This is the sixth representation of the default set, following Bowers
and Burton, "A More Powerful Test for Randomized Experiments" (working
paper), and it reproduces that paper's computation.
