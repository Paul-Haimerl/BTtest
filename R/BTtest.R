#' Barigozzi & Trapani (2022)  Test
#'
#' @description Runs the testing routine proposed in Barigozzi & Trapani (2022) to estimate the number and types of common trends in a nonstationary panel.
#' The method can identify the existence of a common factor subject to a linear trend, as well as the number of zero-mean \eqn{I(1)} and zero-mean \eqn{I(0)} factors.
#'
#' @param X a \eqn{T \times N} numerical \code{matrix} or \code{data.frame} of observations.
#' @param r_max the maximum number of factors to consider. Default is 10. Note that chaning \code{r_max} does not alter the test result for any individual \code{r}.
#' @param alpha the significance level. Default is 0.05.
#' @param BT1 logical. If \code{TRUE}, a less conservative eigenvalue rescaling scheme is used. In small samples, \code{BT1 = FALSE} will result in fewer estimated factors. Default is \code{TRUE}.
#'
#' @details For details on the testing procedure I refer to Barigozzi & Trapani (2022, sec. 4).
#'
#' @examples
#' # Simulate a nonstationary panel
#' X <- sim_DGP(N = 100, n_Periods = 200)
#'
#' # Obtain the estimated number of factors (i) with a linear trend (r_1), (ii) zero-mean I(1) (r_2)
#' # and (iii) zero-mean I(0) (r_3)
#' BTtest(X = X, r_max = 10, alpha = 0.05, BT1 = TRUE)
#' @references Barigozzi, M., & Trapani, L. (2022). Testing for common trends in nonstationary large datasets. *Journal of Business & Economic Statistics*, 40(3), 1107-1122. \doi{10.1080/07350015.2021.1901719}
#'
#' @author Paul Haimerl
#'
#' @return A vector with the estimated number of (i) factors with a linear trend (\eqn{r_1}), (ii) zero-mean \eqn{I(1)} factors (\eqn{r_2}) and (ii) zero-mean \eqn{I(0)} factors (\eqn{r_3}).
#'
#' @export
BTtest <- function(X, r_max = 10, alpha = 0.05, BT1 = TRUE){
  X <- as.matrix(X)
  r_max <- round(r_max)
  if (!is.numeric(X)) stop("`X` must be a matrix or data.frame of numerical values\n")
  if (!is.logical(BT1)) stop("`BT1` must take a logical value\n")
  if (r_max < 1) stop("`r_max` must be a positive integer\n")
  if (alpha <= 0 | alpha >= 1) stop("`alpha` must be a numeric value between 0 and 1\n")
  BTtestRoutine(X = X, r_max = r_max, alpha = alpha, BT1 = BT1)
}


#' Bai (2004) IPC
#'
#' @description Calculates the Integrated Panel Criteria (\emph{IPC}) to estimate the total number of common trends in a nonstationary panel as proposed by Bai (2004).
#'
#' @param X a \eqn{T \times N} numerical \code{matrix} or \code{data.frame} of observations.
#' @param r_max the maximum number of factors to consider. Default is 10.
#'
#' @details For further details on the three criteria and their respective differences, I refer to Bai (2004, sec. 3).
#'
#' @examples
#' # Simulate a nonstationary panel
#' X <- sim_DGP(N = 100, n_Periods = 200)
#'
#' # Obtain the estimated number of common factors pre criterion
#' BaiIPC(X = X, r_max = 10)
#' @references Bai, J. (2004). Estimating cross-section common stochastic trends in nonstationary panel data. *Journal of Econometrics*, 122(1), 137-183. \doi{10.1016/j.jeconom.2003.10.022}
#'
#' @author Paul Haimerl
#'
#' @return A vector of the estimated number of factors for each of the three criteria.
#'
#' @export
BaiIPC <- function(X, r_max = 10){
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("`X` must be a matrix or data.frame of numerical values\n")
  r_max <- round(r_max)
  if (r_max < 1) stop("`r_max` must be a positive integer\n")
  BaiIPCRoutine(X = X, r_max = r_max)
}
