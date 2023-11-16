#' Runs the Barigozzi Trapani test for the number of factors
#'
#' @param X (T x N) matrix of observations
#' @param rmax maximum number of factors to consider
#' @param alpha significance level
#'
#' @return vector with the estimated number of factors
#'
#' @export

BT.test <- function(X, rmax = 20, alpha = .05) return(nFactors(X, rmax, alpha))
