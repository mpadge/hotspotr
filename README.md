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
    ## -------------------------------------|-------------------------------|
    ##  (0.1, 0.1)  0.00    0.29    0.0005  0.7492  |   (-0.01, 0.13)   10  |
    ##  (0.1, 0)    0.00    0.20    0.0001  0.6755  |   (0.06, 0.14)    10  |
    ## ----------------------------------------------------------------------

`run_tests` can also be used to test a neutral field against a statistial ensemble of neutral fields:

``` r
run_tests (size=10, ntests=100, neutral=TRUE)
```

    ##              differences     p-values        |                       |
    ##   alpha      raw     AC      raw     AC      |   alpha           n   |
    ## ---------------------------------------------|-----------------------|
    ##  (0.1, 0.1)  0.00    0.19    0.0329  0.2031  |   (0.55, 0.08)    10  |
    ##  (0.1, 0)    0.00    0.09    0.0263  0.9318  |   (0.50, 0.09)    10  |
    ## ----------------------------------------------------------------------
