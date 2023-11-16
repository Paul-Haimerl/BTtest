
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BTtest

<!-- badges: start -->
<!-- [![CRAN\_Version\_Badge](http://www.r-pkg.org/badges/version/BTtest)](https://cran.r-project.org/package=BTtest) -->
<!-- [![CRAN\_Downloads\_Badge](https://cranlogs.r-pkg.org/badges/grand-total/BTtest)](https://cran.r-project.org/package=BTtest) -->
<!-- badges: end -->

Conveniently execute the testing procedure proposed in Barigozzi &
Trapani ([2022](https://doi.org/10.1080/07350015.2021.1901719)).

## Installation

You can install the development version of BTtest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Paul-Haimerl/BTtest")
```

After installation, attach the package as usual:

``` r
library(BTtest)
```

## Data

The `BTtest` packages includes a function to automatically simulate a
panel subject to common nonstationary trends. Calling this function
provides us with a quick dataset to apply the testing procedure to.

``` r
# Simulate a DGP containing a factor with a linear drift (r1 d1 = 1) and I(1) process (d2 = 1), 
# one zero-mean I(1) factor (r2 = 1) and two zero-mean I(0) factors (r3 = 2)
X <- sim.DGP(N = 100, n.Periods = 200, drift = TRUE, driftI1 = TRUE, r_I1 = 2, r_I0 = 2)
```

For specifics on the DGP I refer to Barigozzi & Trapani
([2022](https://doi.org/10.1080/07350015.2021.1901719), sec. 5).

## The Barigozzi & Trapani test

In order to run the test, one only needs to specify a suitable upper
limit on the number of factors to be tested for as well as a
significance level:

``` r
BT.test_R(X = X, rmax = 10, alpha = 0.05))
```

## References

- Matteo Barigozzi & Lorenzo Trapani (2022) Testing for Common Trends in
  Nonstationary Large Datasets, Journal of Business & Economic
  Statistics, 40:3, 1107-1122, DOI:
  [10.1080/07350015.2021.1901719](https://doi.org/10.1080/07350015.2021.1901719)
