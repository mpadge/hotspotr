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
seed <- 18
ydat <- ives2D (size=10, seed=seed)
x11 ()
par (mfrow=c(1,2))
test <- test1d (ydat, plot=TRUE)
test <- test2d (ydat, plot=TRUE)
```

![](README_files/figure-markdown_github/demo-1.png)

Then the text output of `run_tests` for a random seed giving more typically low p-values

``` r
run_tests (ntests=1000)
```

    ##   dim    |   alpha   diff    p(T)    |   alpha       n   |
    ## --------|-------------------------------|-------------------------------|
    ##   1  |   (0.1, 0.1)  2.13    0.0000  |   (0.11, 0.16)    13  |
    ##   1  |   (0.1, 0)    2.25    0.0000  |   (0.16, 0.00)    10  |
    ##   2  |   (0.1, 0.1)  2.31    0.0000  |   (0.10, 0.11)    8   |
    ##   2  |   (0.1, 0)    2.21    0.0000  |   (0.16, 0.02)    5   |
    ## -------------------------------------------------------------------------
