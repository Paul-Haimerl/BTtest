#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec eigenNorm(const arma::colvec &nu, const bool &BT1)
{
  arma::colvec nu_red;
  if (!BT1)
  {
    nu_red = nu.tail(nu.size() - 1);
  }
  else
  {
    nu_red = nu;
  }
  int N = nu_red.size();
  arma::vec nuSum = arma::sort(arma::cumsum(arma::sort(nu_red)), "descend");
  arma::vec scalingVec = (N - arma::regspace(1, N) + 1);
  arma::vec nuBar = nuSum / scalingVec / 4.0;
  arma::vec output;
  if (!BT1)
  {
    output = arma::join_cols(nuBar, arma::vec{nuBar(N - 1)});
  }
  else
  {
    output = nuBar;
  }
  return output;
}


arma::mat standEigVals(const arma::mat &X, int T, const int &N, const double &delta, const bool &BT1)
{
  arma::mat XtX = trans(X) * X;
  arma::mat DeltaX = diff(X, 1);

  arma::mat Sigma_1 = XtX / pow(T, 3);
  arma::mat Sigma_2 = XtX / pow(T, 2);
  arma::mat Sigma_3 = (trans(DeltaX) * DeltaX) / T;

  // Eq. 9
  arma::colvec nu_1 = sort(arma::eig_sym(Sigma_1), "descend");
  // Eq. 10
  arma::colvec nu_2 = sort(arma::eig_sym(Sigma_2), "descend");
  // Eq. 17
  arma::colvec nu_3 = sort(arma::eig_sym(Sigma_3), "descend");

  // Apply scaling
  // Eq. 18
  arma::vec nu_3_bar = eigenNorm(nu_3, BT1);
  // Eq. 20
  arma::colvec phi_1 = exp(nu_1 / (nu_3_bar * pow(N, delta)));
  // Eq. 25
  arma::colvec phi_2 = exp(nu_2 / (nu_3_bar * pow(N, delta) / log(log(T))));
  arma::colvec phi_3 = exp((nu_3 / 4) / (nu_3_bar * pow(N, delta)));

  arma::mat Phi = join_rows(phi_1, phi_2, phi_3);
  return Phi;
}


float randomTest(const arma::mat &Phi, const int &indx, const int &p, const int &R)
{
  float phi = Phi.col(indx)(p);
  // Set the thresholds according to the zeroes of the Hermite polynomial
  arma::vec U = {-2.4, -.75, .75, 2.4};
  arma::vec Varrho = arma::zeros(4);
  // Compute the corresponding statistics
  arma::vec xi, zeta_tilde;
  arma::uvec zeta;
  float zeta_sum;
  for (int i = 0; i < 4; i++)
  {
    // Step A1.1
    arma::vec xi(R, arma::fill::randn);
    // Step A1.2
    zeta = arma::conv_to<arma::uvec>::from(xi * phi <= U(i));
    zeta_tilde = (arma::conv_to<arma::vec>::from(zeta) - .5) / .5;
    // Step A1.3
    zeta_sum = arma::accu(zeta_tilde);
    Varrho(i) = zeta_sum / sqrt(R);
  }
  arma::vec w = {0.05, 0.45, 0.45, 0.05};
  // Step A1.4
  float Theta = arma::accu(w % pow(Varrho, 2));
  float pvalue = R::pchisq(Theta, 1, 0, 0);
  return pvalue;
}


int randomTestWrapper(const arma::mat &Phi, const int &indx, const int &r_max, const float &sig)
{
  int N = Phi.n_rows;
  int R_prime = fmax(floor(N / 3.0), 100);
  int R;
  int r_hat = 0;
  // Iteratively test the different eigenvalues
  for (int r = 0; r < r_max; r++)
  {
    if (r == 0)
    {
      R = 2 * N;
    }
    else
    {
      R = R_prime;
    }
    float pvalue = randomTest(Phi, indx, r, R);
    if (pvalue > sig)
    {
      r_hat = r + 1;
    }
    else
    {
      break;
    }
  }
  return r_hat;
}


//' Barigozzi & Trapani ([2022](https://doi.org/10.1080/07350015.2021.1901719))  test
//'
//' @description Runs the testing routine proposed in Barigozzi & Trapani ([2022](https://doi.org/10.1080/07350015.2021.1901719)) to estimate the number and types of common trends in a nonstationary panel.
//'
//' @param X a (T x N) matrix of observations
//' @param r_max the maximum number of factors to consider
//' @param alpha the significance level
//' @param BT1 logical. If TRUE, a more conservative eigenvalue rescaling scheme is used
//'
//' @details For details on the testing procedure I refer to Barigozzi & Trapani ([2022](https://doi.org/10.1080/07350015.2021.1901719), sec. 4).
//'
//' @examples
//' # Simulate a nonstationary panel
//' # X <- sim_DGP(N = 100, n_Periods = 200)
//'
//' # Obtain the estimated number of factors (i) with a linear trend (r_1), (ii) zero-mean I(1) (r_2)
//' # and (iii) zero-mean I(0) (r_3)
//' # BTtest(X = X, rmax = 10, alpha = 0.05, BT1 = TRUE)
//' @references Barigozzi, M., & Trapani, L. (2022). Testing for common trends in nonstationary large datasets. *Journal of Business & Economic Statistics*, 40(3), 1107-1122. DOI: [10.1080/07350015.2021.1901719](https://doi.org/10.1080/07350015.2021.1901719)
//'
//' @return A vector with the estimated number of (1) factors with a linear trend (2) zero-mean I(1) factors and (3) zero-mean I(0) factors.
// [[Rcpp::export]]
NumericVector BTtest(const arma::mat &X, const int &r_max, const double &alpha, const bool &BT1)
{

  //------------------------------//
  // Preliminaries                //
  //------------------------------//

  int T = X.n_rows;
  int N = X.n_cols;
  double beta = log(N) / log(T);
  double delta;
  if (beta < .5)
  {
    delta = 1e-5;
  }
  else
  {
    delta = 1e-5 + 1 - .5 / beta;
  }
  double sig = alpha / fmin(N, T);

  //------------------------------//
  // Obtain the eigenvalues       //
  //------------------------------//

  arma::mat Phi = standEigVals(X, T, N, delta, BT1);

  //------------------------------//
  // Tests for r_1, r_star and r  //
  //------------------------------//

  arma::vec r_hat = arma::zeros(3);

  for (int i = 0; i < 3; i++)
  {
    r_hat(i) = randomTestWrapper(Phi, i, r_max, sig);
  }

  if (r_hat(0) >= 1)
  {
    r_hat(0) = 1;
  }

  //------------------------------//
  // Individual factor types      //
  //------------------------------//

  int r_1_hat = r_hat(0);
  int r_2_hat = r_hat(1) - r_1_hat;
  int r_3_hat = r_hat(2) - r_hat(1);

  NumericVector output = NumericVector::create(
    Named("r_1_hat") = r_1_hat,
    Named("r_2_hat") = r_2_hat,
    Named("r_3_hat") = r_3_hat);
  return output;
}

float penal_1(const float &N, const float &T)
{
  float penalty = ((N + T) / (N * T)) * log((N * T) / (N + T));
  return penalty;
}
float penal_2(const float &N, const float &T)
{
  float penalty = ((N + T) / (N * T)) * log(pow(fmin(sqrt(T), sqrt(N)), 2));
  return penalty;
}
arma::vec penal_3(const int &N, const int &T, const int &r_max)
{
  arma::vec penalty = ((N + T - arma::regspace(1, r_max)) / (N * T)) * log(N * T);
  return penalty;
}


//' Bai ([2004](https://doi.org/10.1016/j.jeconom.2003.10.022)) IPC
//'
//' @description Calculates the Integrated Panel Criterions to estimate the number of common trends in a nonstationary panel as proposed in Bai (2004).
//'
//' @param X a (T x N) matrix of observations
//' @param r_max the maximum number of factors to consider
//'
//' @details For further details on the criterion, I refer to Bai ([2004](https://doi.org/10.1016/j.jeconom.2003.10.022), sec. 3).
//'
//' @examples
//' # Simulate a nonstationary panel
//' # X <- sim_DGP(N = 100, n_Periods = 200)
//'
//' # Obtain the estimated number of common factors pre criterion
//' # BaiIPC(X = X, r_max = 10)
//' @references Bai, J. (2004). Estimating cross-section common stochastic trends in nonstationary panel data. *Journal of Econometrics*, 122(1), 137-183. DOI: [10.1016/j.jeconom.2003.10.022](https://doi.org/10.1016/j.jeconom.2003.10.022)
//'
//' @return A vector of the estimated number of factors per criterion.
// [[Rcpp::export]]
NumericVector BaiIPC(const arma::mat &X, const int &r_max)
{

  //------------------------------//
  // Preliminaries                //
  //------------------------------//

  int T = X.n_rows;
  int N = X.n_cols;
  float alpha = T / (4 * log(log(T)));
  arma::mat Xt = trans(X);
  arma::mat Sigma = (Xt * X) / pow(T, 2);

  //------------------------------//
  // Obtain the fitness           //
  //------------------------------//

  arma::mat eigenVec, Resid, Lambda_F, Resid_sq, e;
  arma::vec eigenVal;
  arma::vec ssqVec = arma::zeros(r_max);
  arma::eig_sym(eigenVal, eigenVec, Sigma);

  for (int r = 0; r < r_max; r++)
  {
    e = eigenVec.tail_cols(r + 1);
    Lambda_F = e * trans(e) * Xt;
    // Eq. 4
    Resid = X - trans(Lambda_F);
    Resid_sq = arma::square(Resid);
    ssqVec(r) = arma::mean(arma::mean(Resid_sq));
  }

  //------------------------------//
  // Calculate the ICs            //
  //------------------------------//

  // Scaling parameter
  float sigma_sq = ssqVec(r_max - 1);
  // Compute the constant part of the penalty
  arma::vec penal_fix = arma::regspace(1, r_max) * sigma_sq * alpha;

  // Eq. 12
  arma::vec IPC_1_seq = ssqVec + penal_fix * penal_1(N, T);
  arma::vec IPC_2_seq = ssqVec + penal_fix * penal_2(N, T);
  arma::vec IPC_3_seq = ssqVec + penal_fix % penal_3(N, T, r_max);
  int IPC_1 = arma::index_min(IPC_1_seq) + 1;
  int IPC_2 = arma::index_min(IPC_2_seq) + 1;
  int IPC_3 = arma::index_min(IPC_3_seq) + 1;

  NumericVector output = NumericVector::create(
    Named("IPC_1") = IPC_1,
    Named("IPC_2") = IPC_2,
    Named("IPC_3") = IPC_3);
  return output;
}
