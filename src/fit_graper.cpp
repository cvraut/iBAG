//C Functions calling the main function from classes defined in the included files

#include <RcppArmadillo.h>
#include "graper_linear_dense_ff.h"
#include "graper_linear_dense_nf.h"
#include "graper_linear_sparse_ff_nogamma.h"
#include "graper_linear_sparse_ff.h"
#include "graper_logistic_dense_ff.h"
#include "graper_logistic_dense_nf.h"
#include "graper_logistic_sparse_ff.h"


// [[Rcpp::depends(RcppArmadillo)]]
//Fitting a normal prior model with only partially factorized variational distribution
//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
Rcpp::List graperCpp_dense_nf(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau, double r_tau,
   double d_gamma, double r_gamma, int max_iter, double th, bool calcELB, bool verbose, int freqELB){

    graper MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB);
    List result =MyModel.fitModel();

    return(result); // transforms an arbitrary object into a SEXP.
}


//Fitting a normal prior model with fully factorized variational distribution
//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
Rcpp::List graperCpp_dense_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau, double r_tau,
                            double d_gamma, double r_gamma, int max_iter, double th, bool calcELB, bool verbose, int freqELB, arma::vec mu_init){

  graper_dense_ff MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB, mu_init);
  List result =MyModel.fitModel();

  return(result);
}


//Fitting a spike and slab prior model with fully factorized variational distribution
//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
Rcpp::List graperCpp_sparse_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau, double r_tau,
                          double d_gamma, double r_gamma, double r_pi, double d_pi, int max_iter, double th, bool calcELB, bool verbose,
                             int freqELB, arma::vec mu_init, arma::vec psi_init){

  graper_sparse_ff MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB, verbose, freqELB, mu_init, psi_init);
  List result =MyModel.fitModel();

  return(result);
}

//Fitting a spike and slab prior model with fully factorized variational distribution but no different gamma's per groups
//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
Rcpp::List graperCpp_sparse_ff_nogamma(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup, double d_tau, double r_tau,
                             double d_gamma, double r_gamma, double r_pi, double d_pi, int max_iter, double th, bool calcELB, bool verbose,
                             int freqELB, arma::vec mu_init, arma::vec psi_init){
    
    graper_sparse_ff_nogamma MyModel(X,y,annot,g,NoPerGroup, d_tau, r_tau, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB, verbose, freqELB, mu_init, psi_init);
    List result =MyModel.fitModel();
    
    return(result);
}

//Fitting a normal prior logistic model with only partially factorized variational distribution
//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
Rcpp::List graperCpp_logistic_nf(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup,
                            double d_gamma, double r_gamma, int max_iter, double th,
                            bool calcELB, bool verbose, int freqELB){

  graper_logistic_nf MyModel(X,y,annot,g,NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB);
  List result =MyModel.fitModel();

  return(result);
}


//Fitting a normal prior logistic model with fully factorized variational distribution
// [[Rcpp::export]]
Rcpp::List graperCpp_logistic_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup,
                            double d_gamma, double r_gamma, int max_iter, double th,
                            bool calcELB, bool verbose, int freqELB, arma::vec mu_init, bool intercept ){

  graper_logistic_ff MyModel(X,y,annot,g,NoPerGroup, d_gamma, r_gamma, max_iter, th, calcELB, verbose, freqELB, mu_init, intercept);
  List result =MyModel.fitModel();

  return(result);
}

//Fitting a spike and slab prior logistic model with fully factorized variational distribution
//' @importClassesFrom Matrix dgTMatrix
// [[Rcpp::export]]
Rcpp::List graperCpp_sparse_logistic_ff(arma::mat X, arma::vec y, arma::Row<int> annot, int g, arma::vec NoPerGroup,
                                double d_gamma, double r_gamma, double r_pi, double d_pi, int max_iter, double th,
                                bool calcELB, bool verbose, int freqELB, arma::vec mu_init, arma::vec psi_init, bool intercept){
    
    graper_sparse_logistic_ff MyModel(X,y,annot,g,NoPerGroup, d_gamma, r_gamma, r_pi, d_pi, max_iter, th, calcELB, verbose, freqELB, mu_init, psi_init,intercept);
    List result =MyModel.fitModel();
    
    return(result);
}
