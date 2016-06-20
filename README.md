[![Build Status](https://travis-ci.org/mpadge/hotspotr.svg?branch=master)](https://travis-ci.org/mpadge/hotspotr) [![codecov](https://codecov.io/gh/mpadge/hotspotr/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/hotspotr)

R package for estimating whether the statistical properties of a spatial pattern of hotspots may be reproduced with a simple neutral model. Currently only works for gridded data.

Install
-------

``` r
devtools::install_github ('mpadge/hotspotr')
```

Test
----

First a demonstration with a seed that produces a 2D field very similar to neutral fields.

``` r
plot.new ()
seed <- 18
ymat <- ives2d (size=10, seed=seed)
test <- test2d (ymat, plot=TRUE)
```

![](fig/demo-moran.png)

The default spatial autocorrelation statistic is Moran's I, with results for other statistics in this case notably different. Geary's C:

``` r
plot.new ()
test <- test2d (ymat, plot=TRUE, ac_type='geary')
```

![](fig/demo-geary.png)

And Getis-Ord:

``` r
plot.new ()
test <- test2d (ymat, plot=TRUE, ac_type='getis')
```

![](fig/demo-getis.png)

Then the text output of `run_tests` for a random seed giving more typically low p-values.

``` r
run_tests (size=10, ntests=100)
```

    ##              differences     p-values        |                       |
    ##   alpha      raw     AC      raw     AC      |   alpha           n   |
    ## ---------------------------------------------|-----------------------|
    ##  (0.1, 0.1)  0.21    0.06    0.1267  0.9009  |   (0.81, 0.09)    10  |
    ##  (0.1, 0)    0.18    0.08    0.1623  0.9859  |   (0.86, 0.09)    9   |
    ## ----------------------------------------------------------------------

Note that the model field tested here has a complex *temporal* structure yet a generic---that is, neutral---\*spatial structure, and that the `p-values` reflect these differences.

`run_tests` can also be used to test a neutral field against a statistial ensemble of neutral fields:

``` r
run_tests (size=10, ntests=100, neutral=TRUE)
```

    ##              differences     p-values        |                       |
    ##   alpha      raw     AC      raw     AC      |   alpha           n   |
    ## ---------------------------------------------|-----------------------|
    ##  (0.1, 0.1)  0.07    0.18    0.9041  0.9848  |   (0.59, 0.10)    11  |
    ##  (0.1, 0)    0.08    0.21    0.7236  0.6503  |   (0.14, 0.07)    10  |
    ## ----------------------------------------------------------------------
