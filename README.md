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
    ##   1  |   (0.1, 0.1)  0.47    0.0023  0.0087  |   (0.11, 0.11)    94  |
    ##   1  |   (0, 0.1)    0.34    0.0000  0.0054  |   (0.02, 0.13)    76  |
    ##   2  |   (0.1, 0.1)  1.83    0.0000  0.0000  |   (0.12, 0.09)    124 |
    ##   2  |   (0, 0.1)    1.76    0.0000  0.0002  |   (-0.01, 0.10)   76  |
    ##   2  |   (0.1, 0)    0.17    0.4440  0.7712  |   (0.10, 0.00)    92  |
    ## ---------------------------------------------------------------------------------
