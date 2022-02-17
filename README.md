
# GFLlogit

<!-- badges: start -->
<!-- badges: end -->

The goal of GFLlogit is to ...

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
