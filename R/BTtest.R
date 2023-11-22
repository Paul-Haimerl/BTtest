#' Barigozzi & Trapani ([2022](https://doi.org/10.1080/07350015.2021.1901719))  Test
#'
#' @description Runs the testing routine proposed in Barigozzi & Trapani ([2022](https://doi.org/10.1080/07350015.2021.1901719)) to estimate the number and types of common trends in a nonstationary panel.
#'
#' @param X a (\emph{T x N}) matrix of observations.
#' @param r_max the maximum number of factors to consider.
#' @param alpha the significance level.
#' @param BT1 logical. If \code{TRUE}, a more conservative eigenvalue rescaling scheme is used. Default is \code{TRUE}.
#'
#' @details For details on the testing procedure I refer to Barigozzi & Trapani ([2022](https://doi.org/10.1080/07350015.2021.1901719), sec. 4).
#'
#' @examples
#' # Simulate a nonstationary panel
#' X <- sim_DGP(N = 100, n_Periods = 200)
#'
#' # Obtain the estimated number of factors (i) with a linear trend (r_1), (ii) zero-mean I(1) (r_2)
#' # and (iii) zero-mean I(0) (r_3)
#' BTtest(X = X, r_max = 10, alpha = 0.05, BT1 = TRUE)
#' @references Barigozzi, M., & Trapani, L. (2022). Testing for common trends in nonstationary large datasets. *Journal of Business & Economic Statistics*, 40(3), 1107-1122. DOI: \doi{10.1080/07350015.2021.1901719}
#'
#' @return A vector with the estimated number of (i) factors with a linear trend (r_1), (ii) zero-mean \emph{I(1)} factors  (r_2) and (3) zero-mean \emph{I(0)} factors (r_3).
#'
#' @export
BTtest <- function(X, r_max = 10, alpha = 0.05, BT1 = TRUE){
  BTtestRoutine(X = X, r_max = r_max, alpha = alpha, BT1 = BT1)
}


#' Bai ([2004](https://doi.org/10.1016/j.jeconom.2003.10.022)) IPC
#'
#' @description Calculates the Integrated Panel Criterions (IPC) to estimate the number of common trends in a nonstationary panel as proposed in Bai ([2004](https://doi.org/10.1016/j.jeconom.2003.10.022)).
#'
#' @param X a (T x N) matrix of observations.
#' @param r_max the maximum number of factors to consider.
#'
#' @details For further details on the criterion, I refer to Bai ([2004](https://doi.org/10.1016/j.jeconom.2003.10.022), sec. 3).
#'
#' @examples
#' # Simulate a nonstationary panel
#' X <- sim_DGP(N = 100, n_Periods = 200)
#'
#' # Obtain the estimated number of common factors pre criterion
#' BaiIPC(X = X, r_max = 10)
#' @references Bai, J. (2004). Estimating cross-section common stochastic trends in nonstationary panel data. *Journal of Econometrics*, 122(1), 137-183. DOI: \doi{10.1016/j.jeconom.2003.10.022}
#'
#' @return A vector of the estimated number of factors per criterion.
#'
#' @export
BaiIPC <- function(X, r_max = 10){
  BaiIPCRoutine(X = X, r_max = r_max)
}
