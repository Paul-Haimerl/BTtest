#' Runs the Barigozzi & Trapani test for the number of factors
#'
#' @description Runs the testing routine proposed in Barigozzi & Trapani (2022) to estimate the number of common trends in a nonstationary panel.
#'
#' @param X a (T x N) matrix of observations
#' @param rmax the maximum number of factors to consider
#' @param alpha the significance level
#'
#' @details For details on the testing procedure I refer to Barigozzi & Trapani (2022, sec. 4).
#'
#' @seealso [sim.DGP()]
#'
#' @examples
#' # Simulate some data
#' X <- sim.DGP(N = 100, n.Periods = 200)
#' BT.test_R(X = X, rmax = 10, alpha = 0.05)
#' @references Matteo Barigozzi & Lorenzo Trapani (2022) Testing for Common Trends in Nonstationary Large Datasets, Journal of Businss & Economic Statistics, 40:3, 1107-1122, DOI: 10.1080/07350015.2021.1901719
#'
#' @return A vector with the estimated number of (1) factors with a linear trend (2) zero-mean I(1) factors and (3) zero-mean I(0) factors.
#'
#' @export

BT.test_R <- function(X, rmax = 10, alpha = .05) {
  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  Tt <- NROW(X)
  N <- NCOL(X)
  beta <- log(N) / log(Tt)
  delta <- ifelse(beta < .5, 1e-5, 1e-5 + 1 - .5 / beta)
  XtX <- t(X) %*% X
  DeltaX <- diff(X)
  sig <- alpha / min(N, Tt)


  #------------------------------#
  #### Obtain eigenvalues     ####
  #------------------------------#

  Sigma_1 <- XtX / Tt^3
  Sigma_2 <- XtX / Tt^2
  Sigma_3 <- t(DeltaX) %*% DeltaX / Tt

  # Eq. 9
  nu_1 <- sort(eigen(Sigma_1)$values, decreasing = T)
  # Eq. 10
  nu_2 <- sort(eigen(Sigma_2)$values, decreasing = T)
  # Eq. 17
  nu_3 <- sort(eigen(Sigma_3)$values, decreasing = T)

  # Apply scaling
  # Eq. 18
  nu_3_bar <- eigenNorm_R(nu = nu_3, k = 0)
  # Eq. 20
  phi_1 <- exp(N^(-delta) * nu_1 / nu_3_bar)
  # Eq. 25
  phi_2 <- exp(N^(-delta) * log(log(Tt)) * nu_2 / nu_3_bar)
  phi_3 <- exp(N^(-delta) * nu_3 / nu_3_bar)

  #------------------------------#
  #### Test for r_1           ####
  #------------------------------#

  # Set the threshold to define the simulated Bernoulli sequence (see step A1.2 in sec. 4)
  u <- sqrt(2)
  pvalue_1 <- randomTest_R(phi = phi_1, R = N, p = 1, u = u)
  r_1_hat <- ifelse(pvalue_1 > sig, 1, 0)

  #------------------------------#
  #### Test for r_star        ####
  #------------------------------#

  r_star_hat <- randomTestWrapper_R(phi = phi_2, rmax = rmax, sig = sig, R = N, u = u)
  r_2_hat <- r_star_hat - r_1_hat

  #------------------------------#
  #### Test for r             ####
  #------------------------------#

  r_hat <- randomTestWrapper_R(phi = phi_3, rmax = rmax, sig = sig, R = N, u = u)
  r_3_hat <- r_hat - r_star_hat

  return(c(r_1_hat = r_1_hat, r_2_hat = r_2_hat, r_3_hat = r_3_hat))
}


#' Runs the iterative randomized testing routine
#'
#' @param phi vector of scaled eigenvalues
#' @param rmax maximum number of ordered eigenvalues to consider
#' @param sig significance level
#' @param R Number of simulated random variables per iteration
#' @param u threshold to define a Bernoulli sequence
#'
#' @return estimated number of factors
#'

randomTestWrapper_R <- function(phi, rmax, sig, R, u) {
  r_hat <- 0
  # Iteratively test the different eigenvalues
  for (r in 1:rmax) {
    pvalue <- randomTest_R(phi = phi, p = r, R = R, u = u)
    if (pvalue > sig) {
      r_hat <- r
    } else {
      break
    }
  }
  return(r_hat)
}


#' Runs the randomized testing routine
#'
#' @param phi vector of scaled eigenvalues
#' @param p pth largest eigenvalue to consider
#' @param R Number of simulated random variables
#' @param u threshold to define a Bernoulli sequence
#'
#' @return p-value of the test

randomTest_R <- function(phi, p, R, u) {
  # set.seed(2)
  phi <- phi[p]
  # Step A1.1
  xi <- stats::rnorm(R)
  # Step A1.2
  zeta <- as.numeric(xi * phi <= u)
  # Step A1.3
  varrho <- 1 / sqrt(R) * sum(zeta - .5)
  # Step A1.4
  tstat <- varrho^2
  pvalue <- stats::pchisq(tstat, df = 1, lower.tail = FALSE)
  return(pvalue)
}


#' Obtain a standardization factor for the eigenvalues according to Eq. 18
#'
#' @param nu ordered (descending) vector of eigenvalues
#' @param k censoring threshold. Largest k-1 eigenvalues are omitted
#'
#' @return Scaling factor

eigenNorm_R <- function(nu, k) {
  if (k != 0) nu <- nu[-(1:(k - 1))]
  nu_bar <- 1 / (4 * (length(nu) - k + 1)) * sum(nu)
  return(nu_bar)
}
