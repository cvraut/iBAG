#include <armadillo>
#include  <boost/math/special_functions/digamma.hpp>

# include <iostream>
# include <chrono>
using get_time = std::chrono::steady_clock ;

using  namespace  arma;
using namespace Rcpp;


class  graper_logistic_ff {

private:
  // values remaining constant
  mat  X;
  vec y;
  Row<int> annot;
  int p,n;
  vec NoPerGroup;
  double d_gamma;
  int max_iter;
  double th;
  bool calcELB, verbose, intercept;
  int freqELB;
  List ListXrowSquared;


  // changing values
  double ELB, beta0, EW_beta0, cond_var_beta0, var_beta0;
  sp_mat D_Lambda_xi;
  mat D_Lambda_xi_hat;
  vec cov_beta0;
  vec  alpha_gamma, beta_gamma;
  vec xi, lambda_xi;
  vec XTxi;
  vec Xtyhat, yhat;
  sp_mat Sigma_beta;
  vec sigma2_beta;
  vec mu_beta;
  vec EW_gamma;
  double diff;
  int n_iter;
  vec EW_betasq;
  vec ELB_trace;
public:

  //initaliser list
  graper_logistic_ff(mat X, vec y, Row<int> annot, int g, vec NoPerGroup,
       double d_gamma, double r_gamma, int max_iter, double th, bool calcELB, bool verbose,
       int freqELB, vec mu_init, bool intercept):
  X(X)                                // design matrix
  , y(y)                                // response vector
  , annot(annot)                        // assignement of each feautre to a group
  , p(X.n_cols)                         //number of samples
  , n(X.n_rows)                         //number of samples
  //, g(g)                                // number of groups
  , NoPerGroup(NoPerGroup)               //number of features per group
  , d_gamma(d_gamma)                    // hyperparameters of gamma distribution for gamma
  //, r_gamma(r_gamma)                    // hyperparameters of gamma distribution for gamma
  , max_iter(max_iter)                  // maximal number of iterations
  , th(th)                              //threshold for ELBO to stop iterations
  , calcELB(calcELB)                    //whether to calculate ELBO
  , verbose(verbose)                    //whether to print intermediate messages
  , intercept(intercept)                //whether to use an itercept term in the model
  , freqELB(freqELB)                    // freuqency of ELB calculation: each freqELB-th iteration ELBO is calculated
  , ListXrowSquared(n)
    , ELB(-std::numeric_limits<double>::infinity())                           //evidence lower bound
  , beta0(0)
  , EW_beta0(0)
 , cond_var_beta0(0)
   , var_beta0(0)
 , D_Lambda_xi(speye(n,n))
  , cov_beta0(p)
 , alpha_gamma(g)                      //parameter of gamma distribution for tau (stays same in each iteration)
  , beta_gamma(g)
  , xi(n)                               // variational parameter, initialised to 0, should better be close to yX\beta
  , lambda_xi(n)                       //lambda(xi)
    , Xtyhat(p)
      , yhat(y-0.5)
  , Sigma_beta(speye(p,p))              //diagonal covariance matrix
    , sigma2_beta(p)                      //variance parameter of normal distribution for beta
  , mu_beta(mu_init)                          //initialise by 0
  , EW_gamma(g)                         //initialise by expected value of a gamma distribution, one value per group
  , diff(th+1)                          // to ensure it is larger than th at the beginning
  , n_iter(0)                           // counter of iterations
  , ELB_trace(max_iter)
  { Xtyhat=trans(X)*yhat;
    EW_gamma.fill(r_gamma/d_gamma);
    alpha_gamma=r_gamma+NoPerGroup/2;
    xi.fill(1);
      cov_beta0.fill(0);
    //calculate often used quantities
    for(int i = 0; i< n; i++) {
      ListXrowSquared(i)= X.row(i).t()%X.row(i).t();
      lambda_xi(i) = lambda(xi(i));
      }
    D_Lambda_xi.diag() = lambda_xi;
    if(intercept) D_Lambda_xi_hat = D_Lambda_xi - lambda_xi*lambda_xi.t()/(accu(lambda_xi));
      else D_Lambda_xi_hat = D_Lambda_xi;
	XTxi = X.t()* lambda_xi;
  }

  //helper functions
  double sigmoid (double x){
    return(1/(1+exp(-x)));
  }
  double lambda (double x) {
    if(x!=0) return (1/(2*x)*(sigmoid(x)-0.5));
      else return(1/8);
  }

  //main function: fit  model
  List fitModel() {
    while(n_iter<max_iter && (std::abs(diff)>th || std::isinf(diff) || std::isnan(diff))){
      iterate();
    }

    if(diff<th){
      Rcout << "ELB converged" << endl;
      ELB_trace=ELB_trace(span(0,n_iter-1));
    }
    else{
      Rcout << "Maximum numbers of iterations reached - no convergence or ELB not calculated" << endl;
    }
        List results=List::create(Named("EW_beta")=mu_beta, Named("EW_gamma")=EW_gamma,
                                            Named("ELB")=ELB,
                                    Named("alpha_gamma")=alpha_gamma,
                                    Named("beta_gamma")=beta_gamma, Named("Sigma_beta")=Sigma_beta, Named("ELB_trace")=ELB_trace,
                                    Named("intercept")=EW_beta0);



    return(results);
  }

  //function to do one iteration
  void iterate(){
    n_iter=n_iter+1;                          //increasing counter by 1
    if(verbose) Rcout << "iteration " << n_iter << endl;
      
      update_param_beta();
      update_exp_beta();
      
      update_param_gamma();
      update_exp_gamma();
      
      if(intercept) update_beta0_and_yhat();
      update_param_xi();

    //optional: calculate ELB every freqELB-th step
    if(calcELB & (n_iter%freqELB==0)) calculate_ELBO();
    ELB_trace(n_iter-1)=ELB;

  }


  //function to calculate updated parameters for beta variational distirbution
  void update_param_beta(){
    // if(verbose) Rcout << "Updating beta.." << endl;
    // auto start_beta=get_time::now();

    vec gamma_annot(p);
    mat XTxi_hat = X.t()* D_Lambda_xi_hat;
    for(int i = 0; i< p; i++) {
      gamma_annot(i)=EW_gamma(annot(i)-1);      // minus one as annot starts counting at 1 instead of 0
    }

    // vec term1(p);
    // term1.fill(0);
    // for(int k = 0; k< n; k++) {
    //   vec XrowSquared = ListXrowSquared(k);
    //   term1 = term1+lambda_xi(k) * XrowSquared;
    // }
    // sigma2_beta = 1/(gamma_annot+2*term1);
      
    //mat XTxi_hatX = XTxi_hat*X;
    vec XTxi_hatX_diag(p);
    for(int i = 0; i< p; i++) {
        XTxi_hatX_diag(i) = accu(XTxi_hat.row(i).t()%X.col(i));
      }
    sigma2_beta = 1/(gamma_annot+2*XTxi_hatX_diag);
    Sigma_beta.diag() = sigma2_beta;

    // mat XTxi = X.t() * D_Lambda_xi;
    vec vec1 = X*mu_beta;
   for(int i = 0; i< p; i++){
      double old_mu_i = mu_beta(i);
      mu_beta(i) = sigma2_beta(i) * as_scalar(Xtyhat(i) - 2*(accu(XTxi_hat.row(i).t()%vec1) -
        accu(XTxi_hat.row(i).t()%X.col(i))*mu_beta(i)) );
      vec1 = vec1 + (mu_beta(i) - old_mu_i)*X.col(i);
    }
    // for(int i = 0; i< p; i++){
    //   double old_mu_i = mu_beta(i);
    //   mu_beta(i) = sigma2_beta(i) * as_scalar(Xtyhat(i) - 2*(accu(XTxi.row(i).t()%vec1) -
    //     accu(XTxi.row(i).t()%X.col(i))*mu_beta(i)) );
    //   vec1 = vec1 + (mu_beta(i) - old_mu_i)*X.col(i);
    // }

    // auto time_beta = get_time::now() - start_beta;
    // if(verbose) Rcout<<"Time required:"<<std::chrono::duration_cast<std::chrono::milliseconds>(time_beta).count()<<" ms "<<endl;
  }

  //function to calculate updated parameters for tau variational distribution
  void update_param_xi(){
    for(int i = 0; i< n; i++) {
    vec XrowSquared_i = ListXrowSquared(i);
      //xi(i) = sqrt(as_scalar(X.row(i)*(Sigma_beta + mu_beta*mu_beta.t())*X.row(i).t()));
        double term1 =accu(X.row(i).t()%mu_beta);
      xi(i) = sqrt(as_scalar(accu(XrowSquared_i%sigma2_beta) + (term1+EW_beta0)*(term1+EW_beta0) +var_beta0 + 2* accu(X.row(i).t()%cov_beta0)));
      lambda_xi(i) = lambda(xi(i));
    }
    D_Lambda_xi.diag() = lambda_xi;
    XTxi = X.t()* lambda_xi;

  }

    //function to calculate updated intercept and yhat
    void update_beta0_and_yhat(){
        D_Lambda_xi_hat = D_Lambda_xi - lambda_xi*lambda_xi.t()/(accu(lambda_xi));
        yhat = y-0.5-beta0*2*lambda_xi;

        Xtyhat=trans(X)*yhat;
        double lambda_xi_bar = accu(lambda_xi);

        double y_bar = accu(y-0.5);
        beta0 = accu(y-0.5)/(2*lambda_xi_bar);

        cond_var_beta0 = 1/(2*lambda_xi_bar);
        EW_beta0 = cond_var_beta0*as_scalar(y_bar-2*lambda_xi.t()*X*mu_beta);

        var_beta0 = cond_var_beta0*(1+cond_var_beta0*accu(XTxi%XTxi%mu_beta));
        cov_beta0 = -cond_var_beta0 * XTxi%sigma2_beta;

    }
    
  //function to calculate updated parameters for gamma variational distribution
  void update_param_gamma(){

    beta_gamma.fill(d_gamma);

    for(int i = 0; i< p; i++){
      int k = annot[i]-1;                          // minus one as annot stars counting at 1 instead of 0
      beta_gamma[k]=beta_gamma[k]+0.5*EW_betasq[i];
    }
  }

  //function to update expected values of beta
  void update_exp_beta(){
    EW_betasq=square(mu_beta)+sigma2_beta;
  }


  //function to update expected values of gamma
  void update_exp_gamma(){
    EW_gamma=alpha_gamma/beta_gamma;
  }


  //function to calculate ELBO
  void calculate_ELBO(){
    if(n_iter==1) Rcout<<"ELB not implemented"<<endl;
  }


};
//
