
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BTtest

<!-- badges: start -->

[![CRAN_Version_Badge](http://www.r-pkg.org/badges/version/BTtest)](https://cran.r-project.org/package=BTtest)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/grand-total/BTtest)](https://cran.r-project.org/package=BTtest)
[![License:MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit/)
<!-- badges: end -->

You are analyzing a panel data set and want to determine if the
cross-sectional units share a linear trend as well as any $I(1)$ or
$I(0)$ dynamics?

Conveniently test for the number and type of common factors in large
nonstationary panels using the routine by Barigozzi & Trapani
([2022](https://doi.org/10.1080/07350015.2021.1901719)).

## Installation

You can install the development version of BTtest from
[GitHub](https://github.com/) with:

``` r
# install.packages('devtools')
devtools::install_github('Paul-Haimerl/BTtest')
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo Paul-Haimerl/BTtest@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\phaim\AppData\Local\Temp\RtmpG0G6Pf\remotes53d87db5228b\Paul-Haimerl-BTtest-75180c5/DESCRIPTION' ...     checking for file 'C:\Users\phaim\AppData\Local\Temp\RtmpG0G6Pf\remotes53d87db5228b\Paul-Haimerl-BTtest-75180c5/DESCRIPTION' ...   ✔  checking for file 'C:\Users\phaim\AppData\Local\Temp\RtmpG0G6Pf\remotes53d87db5228b\Paul-Haimerl-BTtest-75180c5/DESCRIPTION' (442ms)
#>       ─  preparing 'BTtest':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#> ─  cleaning src
#>       ─  checking for LF line-endings in source and make files and shell scripts (674ms)
#>       ─  checking for empty or unneeded directories
#>       ─  building 'BTtest_0.10.1.tar.gz'
#>      
#> 
library(BTtest)
```

The stable version is available on CRAN:

``` r
install.packages('BTtest')
```

## Data

The `BTtest` packages includes a function that automatically simulates a
panel with common nonstationary trends:

``` r
set.seed(1)
# Simulate a DGP containing a factor with a linear drift (r1 = 1, d1 = 1 -> drift = TRUE) and 
# I(1) process (d2 = 1 -> drift_I1 = TRUE), one zero-mean I(1) factor 
# (r2 = 1 -> r_I1 = 2; since drift_I1 = TRUE) and two zero-mean I(0) factors (r3 = 2 -> r_I0 = 2)
X <- sim_DGP(N = 100, n_Periods = 200, drift = TRUE, drift_I1 = TRUE, r_I1 = 2, r_I0 = 2)
```

For specifics on the DGP, I refer to Barigozzi & Trapani
([2022](https://doi.org/10.1080/07350015.2021.1901719), sec. 5).

## The Barigozzi & Trapani ([2022](https://doi.org/10.1080/07350015.2021.1901719)) test

To run the test, the user only needs to pass a $T \times N$ data matrix
`X` and specify an upper limit on the number of factors (`r_max`), a
significance level (`alpha`) and whether to use a less (`BT1 = TRUE`) or
more conservative (`BT1 = FALSE`) eigenvalue scaling scheme:

``` r
BTresult <- BTtest(X = X, r_max = 10, alpha = 0.05, BT1 = TRUE)
print(BTresult)
#> r_1_hat r_2_hat r_3_hat 
#>       1       1       2
```

Differences between `BT1 = TRUE/ FALSE`, where `BT1 = TRUE` tends to
identify more factors compared to `BT1 = FALSE`, quickly vanish when the
panel includes more than 200 time periods ([Barigozzi & Trapani
2022](https://doi.org/10.1080/07350015.2021.1901719), sec. 5; [Trapani,
2018](https://doi.org/10.1080/01621459.2017.1328359), sec. 3).

`BTtest` returns a vector indicating the existence of (i) a factor
subject to a linear trend ($r_1$), the number of (ii) zero-mean $I(1)$
factors ($r_2$) and the number of (iii) zero-mean $I(0)$ factors
($r_3$). Note that only one factor with a linear trend can be
identified.

## The Bai ([2004](https://doi.org/10.1016/j.jeconom.2003.10.022)) Integrated Information Criterion

An alternative way of estimating the total number of factors in a
nonstationary panel are the Integrated Information Criteria by Bai
([2004](https://doi.org/10.1016/j.jeconom.2003.10.022)). The package
also contains a function to easily evaluate this measure:

``` r
IPCresult <- BaiIPC(X = X, r_max = 10)
print(IPCresult)
#> IPC_1 IPC_2 IPC_3 
#>     2     2     2
```

## References

- Bai, J. (2004). Estimating cross-section common stochastic trends in
  nonstationary panel data. *Journal of Econometrics*, 122(1), 137-183.
  DOI:
  [10.1016/j.jeconom.2003.10.022](https://doi.org/10.1016/j.jeconom.2003.10.022)
- Barigozzi, M., & Trapani, L. (2022). Testing for common trends in
  nonstationary large datasets. *Journal of Business & Economic
  Statistics*, 40(3), 1107-1122. DOI:
  [10.1080/07350015.2021.1901719](https://doi.org/10.1080/07350015.2021.1901719)
- Trapani, L. (2018). A randomized sequential procedure to determine the
  number of factors. *Journal of the American Statistical Association*,
  113(523), 1341-1349. DOI:
  [10.1080/01621459.2017.1328359](https://doi.org/10.1080/01621459.2017.1328359)
