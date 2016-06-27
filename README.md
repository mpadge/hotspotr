[![Build Status](https://travis-ci.org/mpadge/hotspotr.svg?branch=master)](https://travis-ci.org/mpadge/hotspotr) [![codecov](https://codecov.io/gh/mpadge/hotspotr/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/hotspotr)

R package for estimating whether the statistical properties of a spatial pattern of hotspots may be reproduced with a simple neutral model. Currently only works for gridded data.

Install
-------

``` r
devtools::install_github ('mpadge/hotspotr')
```

Demonstration
----

First a demonstration with a seed that produces a 2D field very similar to neutral fields.

``` r
plot.new ()
seed <- 18
dat <- ives (size=10, seed=seed)
test <- test_hotspots (z=dat$z, nbs=dat$nbs, plot=TRUE)
```

![](fig/demo-moran.png)

The default spatial autocorrelation statistic is Moran's I, with results for other statistics in this case notably different. Geary's C:

``` r
plot.new ()
test <- test_hotspots (z=dat$z, nbs=dat$nbs, plot=TRUE, ac_type='geary')
```

![](fig/demo-geary.png)

And Getis-Ord:

``` r
plot.new ()
test <- test_hotspots (z=dat$z, nbs=dat$nbs, plot=TRUE, ac_type='getis')
```

![](fig/demo-getis.png)

Then the text output of `run_tests` for a random seed giving more typically low p-values.

``` r
dat <- ives (size=10)
run_tests (nbs=dat$nbs, z=dat$z, ntests=100)
```

    ##              differences     p-values        |                       |
    ##   alpha      raw     AC      raw     AC      |   alpha           n   |
    ## ---------------------------------------------|-----------------------|
    ##  (0.1, 0.1)  0.20    0.00    0.8104  0.4635  |   (0.87, -0.05)   10  |
    ##  (0.1, 0)    0.19    0.00    0.8073  0.5280  |   (0.28, 0.02)    10  |
    ## ----------------------------------------------------------------------

Note that the model field tested here has a complex *temporal* structure yet a generic---that is, neutral---\*spatial structure, and that the `p-values` reflect these differences.

`run_tests` can also be used to test a neutral field against a statistial ensemble of neutral fields:

``` r
run_tests (neutral=TRUE)
```

    ##              differences     p-values        |                       |
    ##   alpha      raw     AC      raw     AC      |   alpha           n   |
    ## ---------------------------------------------|-----------------------|
    ##  (0.1, 0.1)  0.29    0.00    0.4908  0.2884  |   (1.13, -0.04)   10  |
    ##  (0.1, 0)    0.28    0.00    0.6964  0.3564  |   (0.38, 0.01)    11  |
    ## ----------------------------------------------------------------------
