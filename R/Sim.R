#' Simulate a DGP according to sec. 5
#'
#' @param N number of cross-sectional units
#' @param Tt number of simulated time periods
#' @param drift logical. If TRUE, a linear trend is included (corresponds to both d_1 and r_1)
#' @param driftI1 logical. If TRUE, an I(1) factor moves around the linear trend. Else an I(0) factor (corresponds to d_2)
#' @param r_I1 total number of non zero-mean I(1) factors (corresponds to r_2 + r_1 * d_2)
#' @param r_I0 number of non zero-mean I(0) factors (corresponds to r_3 + r_1 * (1 - d_2))
#' @param factor logical. If TRUE, the factor matrix is returned. Else the simulated observations
#'
#' @return matrix with the simulated observations
#'
#' @export

sim.DGP <- function(N = 100, Tt = 200, drift = TRUE, driftI1 = TRUE, r_I1 = 2, r_I0 = 1, factor = FALSE) {
  set.seed(2)

  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  # For consistency with Barigozzi and Trapani (2021)
  d_1 <- r_1 <- as.numeric(drift)
  # Note that (r_1 == 1, d_1 == 0, d_2 == 1, r_2 == x) is observationally identical with (r_1 == 0, r_2 == x + 1) (same with d_2 == 0 and r_3)
  # As a consequence, I normalize r_1 == d_1 to save one extra (unnecessary) argument for the function
  d_2 <- as.numeric(driftI1)
  if (driftI1 & drift & r_I1 == 0) stop("Either set r_I1 > 0 or driftI1 == FALSE\n")
  if (!driftI1 & drift & r_I0 == 0) stop("Either set r_I1 > 0 or driftI1 == FALSE\n")
  r_2 <- r_I1 - r_1 * d_2
  r_3 <- r_I0 - r_1 * (1 - d_2)
  r <- r_1 + r_2 + r_3

  #------------------------------#
  #### Factor loadings        ####
  #------------------------------#

  Lambda <- simLambda(r = r, N = N)

  #------------------------------#
  #### Factors                ####
  #------------------------------#

  # Generate the I(1) factors according to random walks with serially correlated unit variances
  I1_fmat <- sapply(1:r_I1, function(x, Tt) simRW(Tt), Tt = Tt)
  if (r_I1 == 0) {
    I1_fmat <- NULL
    F2_indx <- NULL
  } else {
    F2_indx <- 1:r_I1
  }

  # Generate the I(0) factors according to some stationary ARMA process with unit variances
  I0_fmat <- sapply(1:r_I0, function(x, Tt) simARMA(Tt), Tt = Tt)
  if (r_I0 == 0) {
    I0_fmat <- NULL
    F3_indx <- NULL
  } else {
    F3_indx <- (1 + r_I1):(r_I1 + r_I0)
  }

  # Generate a linear trend with slope 1
  if (drift) {
    trend <- 1:Tt
    # Attach the trend to one of the factors and order the factor matrix
    if (driftI1) {
      F_1 <- I1_fmat[, 1] + trend
      I1_fmat <- as.matrix(I1_fmat[, -1])
      F2_indx <- F2_indx[-1]
    } else {
      F_1 <- I0_fmat[, 1] + trend
      I0_fmat <- as.matrix(I0_fmat[, -1])
      F2_indx <- F2_indx + 1
      F3_indx <- F3_indx[-1]
    }
  }

  # Rescale and balance the factors
  # If present, use F1 as the baseline, else the first F_2 factor
  if (drift) {
    norm <- sum((as.matrix(diff(F_1)) %*% t(Lambda[, 1]))^2)
    # Rescale F_2
    F_2 <- rescaleFactors(Fmat = I1_fmat, Lambda = Lambda, norm = norm, indx = F2_indx, fd = TRUE)
    # Rescale F_3
    F_3 <- rescaleFactors(Fmat = I0_fmat, Lambda = Lambda, norm = norm, indx = F3_indx, fd = FALSE)
  } else if (!is.null(F2_indx)) {
    F_1 <- NULL
    F_2 <- I1_fmat
    norm <- sum((as.matrix(diff(I1_fmat)) %*% t(Lambda[, F2_indx]))^2)
    # Rescale F_3
    F_3 <- rescaleFactors(Fmat = I0_fmat, Lambda = Lambda, norm = norm, indx = F3_indx, fd = FALSE)
  } else {
    F_1 <- NULL
    F_2 <- NULL
    F_3 <- I0_fmat
  }
  Fmat <- cbind(F_1, F_2, F_3)
  if (factor) return(Fmat)

  #------------------------------#
  #### Cross-sectional errors ####
  #------------------------------#

  U <- simCrossSecErr(Tt = Tt, N = N)

  # Generate the signal to noise ratio
  theta <- setsignal2noise(U = U, Lambda = Lambda, Fmat = Fmat, drift = drift, F2_indx = F2_indx, F3_indx = F3_indx)

  # Compute the final observation
  X <- Fmat %*% t(Lambda) + sqrt(theta) * U
  return(X)
}


#' Calculate signal to noise ratio
#'
#' @param U (T x N) matrix holding cross-sectional errors
#' @param Lambda (N x r) matrix holding loadings
#' @param Fmat (T x r) matrix holding the factors
#' @param drift logical. If TRUE, a linear trend is included (corresponds to both d_1 and r_1)
#' @param F2_indx vector indicating the columns of f2 factors
#' @param F3_indx vector indicating the columns of f3 factors
#'
#' @return scaling factor governing the signal to noise ratio
#'

setsignal2noise <- function(U, Lambda, Fmat, drift, F2_indx, F3_indx) {
  denom <- sum(diff(U)^2)
  if (drift) {
    F1_indx <- 1
  } else {
    F1_indx <- NULL
  }
  DeltaF <- diff(Fmat)
  F1_sum <- as.matrix(Lambda[, F1_indx]) %*% t(DeltaF[, F1_indx])
  F2_sum <- as.matrix(Lambda[, F2_indx]) %*% t(DeltaF[, F2_indx])
  F3_sum <- as.matrix(Lambda[, F3_indx]) %*% t(DeltaF[, F3_indx])
  num <- sum((F1_sum + F2_sum + F3_sum)^2)
  theta <- .5 * num / denom
  return(theta)
}


#' Calibrates the factors among each other
#'
#' @param Fmat (T x r_i) matrix holding the factors for i = 1, 2, 3
#' @param Lambda (N x r) matrix holding loadings
#' @param norm target sum of squares
#' @param indx vector indicating the columns of relevant factors
#' @param fd logical. If TRUE, the first differences of factors are used
#'
#' @return (T x r_i) matrix with calibrated factors
#'

rescaleFactors <- function(Fmat, Lambda, norm, indx, fd) {
  if (is.null(dim(Fmat))) {
    return(NULL)
  }
  if (fd) {
    factorMat_fd <- diff(Fmat)
  } else {
    factorMat_fd <- Fmat
  }
  ssq <- sum((as.matrix(factorMat_fd) %*% t(Lambda[, indx]))^2)
  scaling <- norm / ssq
  Fmat <- Fmat * sqrt(scaling)
  return(Fmat)
}


#' Draws cross-sectional errors according to Eq. 34
#'
#' @param N number of cross-sectional units
#' @param Tt number of simulated time periods
#' @param a parameter governing the serial correlation
#' @param b parameter governing the cross-sectional correlation
#' @param C parameter governing the extent of the cross-sectional correlation
#'
#' @return (T x N) matrix of cross-sectional errors
#'

simCrossSecErr <- function(Tt, N, a = .5, b = .5, C = min(floor(N / 20), 10)) {
  V <- matrix(rnorm(Tt * N), nr = Tt)
  # Specify the cross-sectional correlation
  k_seq <- 1:C
  # Here deviation from eq. 34, by allowing the i index to circle around after reaching N
  V <- V + b * cbind(V[, -k_seq], V[, k_seq])
  # Specify the serial autocorrelation
  U <- apply(V, 2, function(v) {
    stats::filter(c(0, v), filter = a, method = "recursive", )[-1]
  })
  return(U)
}


#' Draws factor loadings in compliance with Ass. 4
#'
#' @param r bumber of factors to simulate
#' @param N bumber of cross-sectional units
#'
#' @return (N x r) matrix of loadings
#'

simLambda <- function(r, N) {
  A <- matrix(stats::rnorm(r * N), nc = r)
  # Perform QR decomposition to obtain orthonormal matrix
  QR <- qr(A)
  Qmat <- qr.Q(QR)
  # Rescale to a standard normal distribution
  Lambda <- (Qmat - mean(Qmat)) / stats::sd(Qmat)
  return(Lambda)
}


#' Simulates zero-mean I(1) factors according to Eq. 32
#'
#' @param Tt number of simulated time periods
#' @param nBurnin number of burn-in periods for the process
#' @param rholimits vector holding upper and lower bounds for the serial correlation parameter
#' @param sd of the process innovations
#' @return vector with the simulated factor
#'

simRW <- function(Tt, nBurnin = 1, rholimits = c(.4, .8), sd = 1) {
  nBurnin <- max(1, nBurnin)
  Tt_tot <- Tt + nBurnin
  # Draw the serial correlation
  rho <- stats::runif(1, rholimits[1], rholimits[2])
  # Draw the innovations
  e_tilde <- stats::rnorm(Tt_tot, sd = sd)
  e <- stats::filter(c(0, e_tilde), filter = rho, method = "recursive", )[-1]
  # Construct the RW
  RW <- cumsum(e)
  return(RW[-(1:nBurnin)])
}


#' Simulates zero-mean I(0) factors according to Eq. 33 with the extension to arbitrary stationary ARMA process
#'
#' @param Tt number of simulated time periods
#' @param pqmax vector holding upper bounds for the AR and MA lag order
#' @param sd of the process innovations
#'
#' @return vector with the simulated factor
#'

simARMA <- function(Tt, pqmax = c(1, 0), sd = 1) {
  # Pick the lag order
  pq <- c(sample(1:pqmax[1], 1), sample(1:pqmax[2], 1))
  # Draw the ar and ma polynomials
  coef <- simARMACoef(pq = pq)
  # Construct the ARMA process
  arma <- as.numeric(stats::arima.sim(list(order = c(pq[1], 0, pq[2]), ar = coef$ar, ma = coef$ma), sd = sd, n = Tt))
  return(arma)
}


#' Simulates the AR and MA coefficients
#'
#' @param pq vector holding the lag orders
#' @param limits vector holding upper and lower bounds for the coefficients
#'
#' @return vector with stationary simulated coefficients
#'

simARMACoef <- function(pq, limits = c(-.5, .5)) {
  # AR
  if (pq[1] == 0) {
    ar <- NULL
  } else {
    while (TRUE) {
      ar <- stats::runif(pq[1], limits[1], limits[2])
      minroots <- min(Mod(polyroot(c(1, -ar))))
      if (minroots > 1) break
    }
    if (!is.null(ar)) names(ar) <- paste0("ar", 1:pq[1])
  }
  # MA
  if (pq[2] == 0) {
    ma <- NULL
  } else {
    ma <- stats::runif(pq[2], -1, 1)
  }
  if (!is.null(ma)) names(ma) <- paste0("ma", 1:pq[2])
  return(list(ar = ar, ma = ma))
}
