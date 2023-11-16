#include <RcppArmadillo.h>
#include <cmath>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


//' Obtain a standardization factor for the eigenvalues according to Eq. 18
//'
//' @param nu ordered (descending) vector of eigenvalues
//' @param k censoring threshold. Largest k-1 eigenvalues are omitted
//'
//' @return Scaling factor

 double EigenNorm(colvec nu, int k){
     colvec nu_red;
     if (k != 0){
         nu_red = nu.tail(nu.size() - k + 1);
     } else {
         nu_red = nu;
     }
     double nu_bar = sum(nu_red) / (4 * (nu_red.size() - k + 1));
     return nu_bar;
 }


//' Obtain a matrix of properly scaled eigenvalues
//'
//' @param X (T x N) matrix of observations
//' @param T number of time periods
//' @param N number of cross-sectional units
//' @param delta scaling factor for N
//'
//' @return (N x r) matrix of scaled eigenvalues

 mat StandEigVals (mat X, int T, int N, double delta){
     mat XtX = trans(X) * X;
     mat DeltaX = diff(X, 1);

     mat Sigma_1 = XtX / pow(T, 3);
     mat Sigma_2 = XtX / pow(T, 2);
     mat Sigma_3 = (trans(DeltaX) * DeltaX) / T;

     // Eq. 9
     colvec nu_1 = sort(eig_sym(Sigma_1), "descend");
     // Eq. 10
     colvec nu_2 = sort(eig_sym(Sigma_2), "descend");
     // Eq. 17
     colvec nu_3 = sort(eig_sym(Sigma_3), "descend");

     // Apply scaling
     // Eq. 18
     double nu_3_bar = EigenNorm(nu_3, 0);
     // Eq. 20
     colvec phi_1 = exp(pow(N, -delta) * nu_1 / nu_3_bar);
     // Eq. 25
     colvec phi_2 = exp(pow(N, -delta) * log(log(T)) * nu_2 / nu_3_bar);
     colvec phi_3 = exp(pow(N, -delta) * nu_3 / nu_3_bar);

     mat Phi = join_rows(phi_1, phi_2, phi_3);
     return Phi;
 }


//' @description Runs the randomized testing routine
//'
//' @param phi vector of scaled eigenvalues
//' @param indx for the factor group in question
//' @param p pth largest eigenvalue to consider
//' @param R Number of simulated random variables
//' @param u threshold to define a Bernoulli sequence
//'
//' @return p-value of the test

 float RandomTest (mat Phi, int indx, int p, int R, float u){
     float phi = Phi.col(indx - 1)(p - 1);
     // Step A1.1
     vec xi(R, fill::randn);
     // Step A1.2
     uvec zeta = conv_to<uvec>::from(xi * phi <= u);
     // Step A1.3
     float zeta_sum = accu(zeta);
     float varrho = 1.0 / sqrt(R) * (zeta_sum - R * .5);
     // Step A1.4
     float tstat = pow(varrho, 2);
     float pvalue = R::pchisq(tstat, 1, 0, 0);
     return pvalue;
 }


//' @description Runs the iterative randomized testing routine
//'
//' @param phi vector of scaled eigenvalues
//' @param rmax maximum number of ordered eigenvalues to consider
//' @param sig significance level
//' @param R Number of simulated random variables per iteration
//' @param u threshold to define a Bernoulli sequence
//'
//' @return estimated number of factors

 int RandomTestWrapper(mat Phi, int indx, int rmax, float sig, int R, float u){
     int r_hat = 0;
     // Iteratively test the different eigenvalues
     for (int r = 1; r <= rmax; r++) {
         float pvalue = RandomTest(Phi, indx, r, R, u);
         if (pvalue > sig) {
             r_hat = r;
         } else {
             break;
         }
     }
     return r_hat;
 }


//' @description Runs the Barigozzi Trapani test for the number of factors
//'
//' @param X (T x N) matrix of observations
//' @param rmax maximum number of factors to consider
//' @param alpha significance level
//'
//' @return vector with the estimated number of factors
// [[Rcpp::export]]

NumericVector nFactors (mat X, int rmax, double alpha){

//------------------------------//
// Preliminaries                //
//------------------------------//

    int T = X.n_rows;
    int N = X.n_cols;
    double beta = log(N) / log(T);
    double delta;
    if (beta < .5){
        delta = 1e-5;
    } else {
        delta = 1e-5 + 1 - .5 / beta;
    }
    double sig = alpha / fmin(N, T);
    double u = sqrt(2);

    //------------------------------//
    // Obtain the eigenvalues       //
    //------------------------------//

    mat Phi = StandEigVals(X, T, N, delta);

    //------------------------------//
    // Test for r_1                 //
    //------------------------------//

    float pvalue_1 = RandomTest(Phi, 1, 1, N, u);
    int r_1_hat;
    if(pvalue_1 > sig){
        r_1_hat = 1;
    } else {
        r_1_hat = 0;
    }

    //------------------------------//
    // Test for r_star              //
    //------------------------------//

    int r_star_hat = RandomTestWrapper(Phi, 2, rmax, sig, N, u);
    int r_2_hat = r_star_hat - r_1_hat;

    //------------------------------//
    // Test for r                   //
    //------------------------------//

    int r_hat = RandomTestWrapper(Phi, 3, rmax, sig, N, u);
    int r_3_hat = r_hat - r_star_hat;

    NumericVector output = NumericVector::create(
       Named("r_1_hat") = r_1_hat,
       Named("r_2_hat") = r_2_hat,
       Named("r_3_hat") = r_3_hat
    );
    return output;
}


