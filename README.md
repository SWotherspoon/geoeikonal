
# geoeikonal

<!-- badges: start -->
<!-- badges: end -->

The `geoeikonal` package computes generalized distances and arrival times on a geographic grid by
solving the eikonal equation by the fast marching algorithm of Sethian (1996).

## Installation

You can install the development version of geoeikonal like so:

```r
# install.packages("remotes")
remotes::install_github("SWotherspoon/geoeikonal")
```


## Example


Reproduce an example from the `terra` package
``` r
library(terra)
library(geoeikonal)
## This is the example for terra::gridDist
r <- rast(ncol=10,nrow=10, vals=1)
r[48] <- 0
r[66:68] <- NA
d <- gridDist(r) 
plot(d)
## With geoeikonal
d <- rsteikonal(r) 
plot(d)
```

