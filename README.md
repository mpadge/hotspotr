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
    ##  (0.1, 0.1)  1.17    0.21    0.0002  0.5559  |   (0.31, 0.17)    10  |
    ##  (0.1, 0)    1.13    0.25    0.0003  0.1035  |   (0.23, 0.11)    10  |
    ## ----------------------------------------------------------------------

Note that the model field tested here has a complex *temporal* structure yet a generic---that is, neutral---\*spatial structure, and that the `p-values` reflect these differences.

`run_tests` can also be used to test a neutral field against a statistial ensemble of neutral fields:

``` r
run_tests (size=10, ntests=100, neutral=TRUE)
```

    ##              differences     p-values        |                       |
    ##   alpha      raw     AC      raw     AC      |   alpha           n   |
    ## ---------------------------------------------|-----------------------|
    ##  (0.1, 0.1)  0.16    0.13    0.1797  0.4289  |   (0.11, 0.12)    10  |
    ##  (0.1, 0)    0.10    0.12    0.3178  0.5556  |   (0.08, 0.11)    10  |
    ## ----------------------------------------------------------------------
