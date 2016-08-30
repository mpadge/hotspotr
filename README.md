[![Build Status](https://travis-ci.org/mpadge/hotspotr.svg?branch=master)](https://travis-ci.org/mpadge/hotspotr) [![codecov](https://codecov.io/gh/mpadge/hotspotr/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/hotspotr)

hotspotr (... work in progress ...)
===================================

There are a number of statistical tools for evaluating the significance of observed spatial patterns of hotspots (such as those in [`spdep`](https://cran.r-project.org/package=spdep)). While such tools enable hotspots to be identified within a given set of spatial data, they do not allow quantification of whether the entire data set in fact reflects significant spatial structuring. For example, numbers of locally significant hotspots must be expected to increase with numbers of spatial observations.

The `R` package `hotspotr` enables the global significance of an observed spatial data set to be quantified by comparing both raw values and their local spatial relationships with those generated from a neutral model. If the global pattern of observed values is able to be reproduced by a neutral model, then any local hotspots may not be presumed significant regardless of the values of local statistics. Conversely, if the global pattern is unable to be reproduced by a neutral model, then local hotspots may indeed be presumed to be statistically significant.

The package is inspired by the work of *Brown, Mehlman, & Stevens* (Ecology 1995) and *Ives & Klopfer* (Ecology 1997). `hotspotr` follows the same premises as these two papers, in examining the extent to which rank--scale distributions can be reproduced by simple neutral models. `hotspotr` compares rank--scale distributions not only of the data of interest, but of corresponding local autocorrelation statistics.

Analysis involves first fitting a model using the function `fit_hotspot_model`, and then testing the significance of that using the function `p-values`.

The remainder of this README documents various exploratory and development phases ...

------------------------------------------------------------------------

Contents
========

[1. Parallel and Rcpp Tests](#1-parallel)

[2. Tests](#2-tests)

------------------------------------------------------------------------

Install
-------

``` r
devtools::install_github ('mpadge/hotspotr')
```

------------------------------------------------------------------------

<a name="1-parallel"></a>1. Parallel and Rcpp Tests
===================================================

Fist set up grid and list of spaital neighbours

``` r
ntests <- 10000
size <- 10
xy <- cbind (rep (seq (size), each=size), rep (seq (size), size))
dhi <- 1 # for rook; dhi=1.5 for queen
nbs <- spdep::dnearneigh (xy, 0, dhi)
```

The function `neutral_hotspots_ntests2` has a `parallel` option that determines whether the tests are conducting using an `Rcpp` (non-parallel) loop, or using `R:parellel`. The latter may be faster under some circumstances for large numbers of tests, although it will definitely be slower for smaller numbers because of the time needed to establish the parallel clusters.

``` r
set.seed (1)
st1 <- system.time (dat1 <- neutral_hotspots (nbs, ntests=ntests))
set.seed (1)
st2 <- system.time (dat2 <- neutral_hotspots (nbs, ntests=ntests,
                                              parallel=TRUE))
st1; st2
```

    ##    user  system elapsed 
    ##   1.396   0.004   1.399

    ##    user  system elapsed 
    ##   0.704   0.044   4.674

The parallel versions do not of course generate identical results, because each core starts with its own random seed, but nevertheless after

``` r
ntests
```

    ## [1] 10000

the differences are very small:

``` r
max (abs (dat1 - dat2))
```

    ## [1] 0.002884018

------------------------------------------------------------------------

<a name="2-tests"></a>2. Tests
==============================

Test the distributional properties of the `meuse` data from the `sp` package, which contain topsoil heavy metal concentrations near Meuse, NL.

``` r
data (meuse, package='sp')
names (meuse)
```

    ##  [1] "x"       "y"       "cadmium" "copper"  "lead"    "zinc"    "elev"   
    ##  [8] "dist"    "om"      "ffreq"   "soil"    "lime"    "landuse" "dist.m"

The function `test_hotspots` requires data to be tested, a list of neighbours (constructed here as `knearneigh`), and a matching list of weights (constructed here as inverse-distance weights):

``` r
xy <- cbind (meuse$x, meuse$y)
nbs <- spdep::knn2nb (spdep::knearneigh (xy, k=4))
dists <- spdep::nbdists (nbs, xy)
# wts calculated in 2 steps to make it explicit
d1 <- lapply (dists, function (i) 1/i)
wts <- lapply (d1, function (i) i / sum (i))
```

Spatial patterns for the different metals can then be statistically compared with neutral models:

``` r
analyse <- function (metal='copper', ntests=1000)
{
    z <- meuse [metal] [,1]
    mod <- fit_hotspot_model (z=z, nbs=nbs, wts=wts, ntests=ntests, verbose=FALSE)
    cat (mod$sd0, ',', mod$ac, ',', mod$niters, '\n')
    p_values (z=z, nbs=nbs, sd0=mod$sd0, alpha=mod$ac, niters=mod$niters,
              ntests=ntests, plot=TRUE)
}
```

For demonstration purposes, `ntests=1000` is sufficient, but larger values will generate more reliable estimates. These functions can be quite time-consuming.

``` r
analyse (metal='cadmium', ntests=1000)
```

    ## stopping search because error is increasing

    ## 0.000001 , 0.2744361 , 6

![](README_files/figure-markdown_github/meuse-cadmium-1.png)![](README_files/figure-markdown_github/meuse-cadmium-2.png)

    ## $p_z
    ## [1] 0
    ## 
    ## $p_ac
    ## [1] 0

![](fig/meuse-cadmium.png)

``` r
analyse (metal='copper', ntests=1000)
```

    ## 0.000001 , 0.593985 , 2

![](README_files/figure-markdown_github/meuse-copper-1.png)![](README_files/figure-markdown_github/meuse-copper-2.png)

    ## $p_z
    ## [1] 0
    ## 
    ## $p_ac
    ## [1] 0.9669162

![](fig/meuse-copper.png)

``` r
analyse (metal='lead', ntests=1000)
```

    ## 0.000001 , 0.2132116 , 6

![](README_files/figure-markdown_github/meuse-lead-1.png)![](README_files/figure-markdown_github/meuse-lead-2.png)

    ## $p_z
    ## [1] 0
    ## 
    ## $p_ac
    ## [1] 0.9999729

![](fig/meuse-lead.png)

``` r
analyse (metal='zinc', ntests=1000)
```

    ## stopping search because error is increasing

    ## 0.000001 , 0.792696 , 1

![](README_files/figure-markdown_github/meuse-zinc-1.png)![](README_files/figure-markdown_github/meuse-zinc-2.png)

    ## $p_z
    ## [1] 0
    ## 
    ## $p_ac
    ## [1] 0.2105138

![](fig/meuse-zinc.png)

These plots demonstrate that in all cases, the observed values themselves (`z` in the figures) can not be reproduced by a neutral model, yet the actual spatial relationships between them can. This indicates that the generating processes can be modelled by assuming nothing other than simple spatial autocorrelation acting on a more complex, non-spatial process.
