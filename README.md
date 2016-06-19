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

Then the text output of `run_tests` for a random seed giving more typically low p-values

``` r
run_tests (size=10, ntests=100)
```

    ##              differences     p-values        |                       |
    ##   alpha      raw     AC      raw     AC      |   alpha           n   |
    ## ---------------------------------------------|-----------------------|
    ##  (0.1, 0.1)  1.16    0.21    0.0003  0.6261  |   (0.10, 0.13)    10  |
    ##  (0.1, 0)    0.99    0.11    0.0006  0.9519  |   (0.38, 0.14)    10  |
    ## ----------------------------------------------------------------------

`run_tests` can also be used to test a neutral field against a statistial ensemble of neutral fields:

``` r
run_tests (size=10, ntests=100, neutral=TRUE)
```

    ##              differences     p-values        |                       |
    ##   alpha      raw     AC      raw     AC      |   alpha           n   |
    ## ---------------------------------------------|-----------------------|
    ##  (0.1, 0.1)  0.34    0.55    0.0530  0.6961  |   (0.34, 0.11)    11  |
    ##  (0.1, 0)    0.48    0.43    0.0222  0.8479  |   (0.94, 0.20)    10  |
    ## ----------------------------------------------------------------------
