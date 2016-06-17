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

![](fig/demo.png)

The default spatial autocorrelation statistic is Moran's I, with results for other statistics in this case notably different. Geary's C:

``` r
plot.new ()
test <- test2d (ymat, plot=TRUE, actype='geary')
```

![](fig/demo-geary.png)

And Getis-Ord:

``` r
plot.new ()
test <- test2d (ymat, plot=TRUE, actype='getis')
```

![](fig/demo-getis.png)

Then the text output of `run_tests` for a random seed giving more typically low p-values

``` r
run_tests (size=10, ntests=1000)
```

    ##      |               differences     p-values        |                       |
    ##  dim |  alpha        raw     AC      raw     AC      |   alpha       n       |
    ## -----|-----------------------------------------------|-----------------------|
    ##   1  |   (0.1, 0.1)  0.00    9.84    0.0000  0.0000  |   (0.12, 0.13)    8   |
    ##   1  |   (0.1, 0)    0.00    9.72    0.0000  0.0000  |   (0.12, 0.00)    10  |
    ##   2  |   (0.1, 0.1)  0.00    10.69   0.0000  0.0000  |   (0.08, 0.11)    13  |
    ##   2  |   (0.1, 0)    0.00    10.39   0.0000  0.0000  |   (0.13, 0.03)    5   |
    ## ------------------------------------------------------------------------------

`run_tests` can also be used to test a neutral field against a statistial ensemble of neutral fields:

``` r
run_tests (size=10, ntests=1000, neutral=TRUE)
```

    ##      |               differences     p-values        |                       |
    ##  dim |  alpha        raw     AC      raw     AC      |   alpha           n   |
    ## -----|-----------------------------------------------|-----------------------|
    ##   1  |   (0.1, 0.1)  0.00    1.46    0.0186  0.0004  |   (0.11, 0.11)    8   |
    ##   1  |   (0.1, 0)    0.00    1.20    0.0283  0.0014  |   (0.10, 0.00)    8   |
    ##   2  |   (0.1, 0.1)  0.00    1.03    0.0190  0.0034  |   (0.11, 0.11)    5   |
    ##   2  |   (0.1, 0)    0.00    0.89    0.0268  0.0095  |   (0.11, -0.01)   5   |
    ## ------------------------------------------------------------------------------
