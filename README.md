[![Build Status](https://travis-ci.org/mpadge/hotspotr.svg?branch=master)](https://travis-ci.org/mpadge/hotspotr) [![codecov](https://codecov.io/gh/mpadge/hotspotr/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/hotspotr)

R package for estimating whether the statistical properties of a spatial pattern of hotspots may be reproduced with a simple neutral model. Currently only works for gridded data.

Install
-------

``` r
devtools::install_github ('mpadge/hotspotr')
```

Test
----

Just a demonstration ....

``` r
run_tests (ntests=10)
```

    ##   dim    |   alpha   diff    p(w)    p(T)    |   alpha       n   |
    ## --------|---------------------------------------|-------------------------------|
    ##   1  |   (0.1, 0.1)  0.11    0.0000  0.4681  |   (0.13, 0.12)    91  |
    ##   1  |   (0.1, 0)    0.05    0.0000  0.6993  |   (0.11, -0.00)   122 |
    ##   2  |   (0.1, 0.1)  0.68    0.5611  0.9908  |   (0.11, 0.11)    124 |
    ##   2  |   (0.1, 0)    0.08    0.0154  0.7657  |   (0.09, 0.01)    75  |
    ## ---------------------------------------------------------------------------------
