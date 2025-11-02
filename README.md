
# fssg: Flexsurv Shotgun Approach

<!-- badges: start -->
<!-- badges: end -->

This package is designed to provide a simple, one-line-of-code approach to testing a variety of parametric survival distributions on a data set. 
This is done by running 60+ different parametric survival curves in order to see which one may be best suited to describing your data! 

## Installation

You can install the development version of fssg like so:

``` r
devtools::install_github('jmrothen/survshotgunR')
```

## Example

The package is designed to be mostly contained to one function, fssg. Simply provide a survival formula, and receive a table of each model run, and how it fits:

``` r
library(survshotgunR)

# sample dataset available in <survival>
library(survival)
surv_shotgun(Surv(time, status)~1, data=aml, dump_models = TRUE)

```


## Developer's Notes

This is still a very young package, so please expect some errors to occur. Please feel free to submit issues via github if you find them.
