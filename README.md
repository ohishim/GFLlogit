
# GFLlogit (v1.0.0)

<!-- badges: start -->
<!-- badges: end -->

**cite this package**:  
Ohishi, M. (2021).
GFLlogit: Coordinate optimization for GFL logistic regression.
R package version 1.0.0.
https://github.com/ohishim/GFLlogit

## Installation

You can install the development version of GFLlogit from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ohishim/GFLlogit")
```

## Example

This package has an example dataset as the two objects `exdata` and `adj`.
The `exdata` has 200 individuals for 10 areas and 20 time points.
The `adj` has adjacent information of 200 spatio-temporal points.
This package estimates spatio-temporal trends.  

The GFLlogit procedure can be executed as follows:

``` r
library(GFLlogit)

m <- exdata$m
y <- exdata$y
adjCD <- split(adj$adj, adj$lab)

res <- GFLlogit(m, y, adjCD)

```

The output `res` takes a form of a list object.
For example, you can obtain GFL estimates by `res$mu.hat`.

## Reference

1. Ohishi, M., Yamamura, M. & Yanagihara, H. (2022).
Coordinate descent algorithm of generalized fused Lasso logistic regression for multivariate trend filtering.
*Jpn. J. Stat. Data Sci.*, **5**, 535-551.
doi: [10.1007/s42081-022-00162-2](https://doi.org/10.1007/s42081-022-00162-2)
