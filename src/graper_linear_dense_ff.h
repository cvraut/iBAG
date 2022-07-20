#include <armadillo>
#include  <boost/math/special_functions/digamma.hpp>

# include <iostream>
# include <chrono>
using get_time = std::chrono::steady_clock ;

using  namespace  arma;
using namespace Rcpp;


class  graper_dense_ff {

private:
  // values remaining constant
  mat  X, XtX;
  vec y, Xty, diagXtX;
  rowvec ytX;
  Row<int> annot;
  double yty;
  int p,n,g;
  vec NoPerGroup;
  double d_tau, r_tau, d_gamma;
  vec r_gamma_mult;
  int max_iter;
  double th;
  bool calcELB, verbose;
  int freqELB;

  // changing values
  double EW_tau, ELB;
  double alpha_tau, beta_tau;
  vec  alpha_gamma, beta_gamma;
  vec sigma2_beta;
  sp_mat Sigma_beta;
  vec mu_beta;
  vec EW_gamma;
  double diff;
  int n_iter;
  vec EW_betasq;
  double EW_leastSquares;
  vec ELB_trace;
public:

  //initaliser list
  graper_dense_ff(mat X, vec y, Row<int> annot, int g, vec NoPerGroup, double d_tau, double r_tau,
       double d_gamma, double r_gamma, int max_iter, double th, bool calcELB, bool verbose,
       int freqELB, vec mu_init):
  X(X)                                // design matrix
  , XtX(trans(X)*X)
  , y(y)                                // response vector
  , Xty(trans(X)*y)
  , diagXtX(XtX.diag())
  , ytX(trans(y)*X)
  , annot(annot)                        // assignement of each feautre to a group
  , yty(as_scalar(trans(y)*y))
  , p(X.n_cols)                         //number of samples
  , n(X.n_rows)                         //number of samples
  , g(g)                                // number of groups
  , NoPerGroup(NoPerGroup)               //number of features per group
  , d_tau(d_tau)                        // hyperparameters of gamma distribution for tau
  , r_tau(r_tau)                        // hyperparameters of gamma distribution for tau
  , d_gamma(d_gamma)                    // hyperparameters of gamma distribution for gamma
//  , r_gamma(r_gamma)                    // hyperparameters of gamma distribution for gamma
  , r_gamma_mult(g)                    // hyperparameters of gamma distribution for gamma
  , max_iter(max_iter)                  // maximal number of iterations
  , th(th)                              //threshold for ELBO to stop iterations
  , calcELB(calcELB)                    //whether to calculate ELBO
  , verbose(verbose)                    //whether to print intermediate messages
  , freqELB(freqELB)                    // freuqency of ELB calculation: each freqELB-th iteration ELBO is calculated
  , EW_tau(r_tau/d_tau)                 //initialise by expected value of a gamma distribution
  , ELB(-std::numeric_limits<double>::infinity())                           //evidence lower bound
  , alpha_tau(r_tau+n/2)                //parameter of gamma distribution for tau (stays same in each iteration)
  , alpha_gamma(g)                      //parameter of gamma distribution for tau (stays same in each iteration)
  , beta_gamma(g)
  , sigma2_beta(p)                      //variance parameter of normal distribution for beta
  , Sigma_beta(speye(p,p))              //diagonal covariance matrix
  , mu_beta(mu_init)                    //initialised randomly
  , EW_gamma(g)                         //initialise by expected value of a gamma distribution, one value per group
  , diff(th+1)                          // to ensure it is larger than th at the beginning
  , n_iter(0)                           // counter of iterations
  , ELB_trace(max_iter)
  {
      vec multgamma = ones<vec>(g); //NoPerGroup;
      r_gamma_mult = r_gamma * multgamma;
      //r_gamma_mult.fill(r_gamma);
      EW_gamma =r_gamma_mult/d_gamma;
      //EW_gamma.fill(r_gamma_mult/d_gamma);
      alpha_gamma=r_gamma_mult+NoPerGroup/2;
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

    List results=List::create(Named("EW_beta")=mu_beta, Named("EW_gamma")=EW_gamma,Named("EW_tau")=EW_tau, Named("ELB")=ELB,
                                    Named("alpha_gamma")=alpha_gamma, Named("alpha_tau")=alpha_tau, Named("beta_tau")=beta_tau,
                                    Named("beta_gamma")=beta_gamma, Named("Sigma_beta")=Sigma_beta, Named("ELB_trace")=ELB_trace);


    return(results);
  }

  //function to do one iteration
  void iterate(){
    n_iter=n_iter+1;                          //increasing counter by 1
    if(verbose) Rcout << "iteration " << n_iter << endl;

      update_param_beta();  
      update_exp_beta();
      update_param_tau();
      update_exp_tau();
      update_param_gamma();
      update_exp_gamma();

    //calculate ELB every freqELB-th step to monitor convergence
    if(calcELB && (n_iter%freqELB==0)) calculate_ELBO();
    ELB_trace(n_iter-1)=ELB;

  }


  //function to calculate updated parameters for beta variational distirbution
  void update_param_beta(){
    vec gamma_annot(p);
    for(int i = 0; i< p; i++) {
      gamma_annot(i)=EW_gamma(annot(i)-1);      // minus one as annot starts counting at 1 instead of 0
    }

      sigma2_beta=1/(EW_tau*diagXtX+gamma_annot);
      Sigma_beta.diag() = sigma2_beta;

    vec vec1 = X*mu_beta;
    for(int i = 0; i< p; i++){
    //mu_beta(i)=sigma2_beta(i)* EW_tau * (Xty(i)- accu(XtX.row(i)%trans(mu_beta)) + XtX(i,i)*mu_beta(i));
    // keep track of old mu for efficient update of vec1
    double old_mu_i = mu_beta(i);
    //update mean mu
    mu_beta(i)=sigma2_beta(i)* EW_tau * (Xty(i)- accu(X.col(i) % vec1) + diagXtX(i)*mu_beta(i));
    //update vec1 (only in new coordinate of mu to avoid recomputing the full pxp product and get linear complexity in inner loops)
    vec1 = vec1 + (mu_beta(i) - old_mu_i)*X.col(i);
    }
  }

  //function to calculate updated parameters for tau variational distribution
  void update_param_tau(){
    beta_tau=d_tau+0.5*EW_leastSquares;
  }

  //function to calculate updated parameters for gamma variational distribution
  void update_param_gamma(){
    beta_gamma.fill(d_gamma);

    for(int i = 0; i< p; i++){
      int k = annot[i]-1;                          // minus one as annot stars counting at 1 instead of 0
      beta_gamma[k]=beta_gamma[k]+0.5*EW_betasq[i];
    }
  }

  //function to update expected values involving beta
  void update_exp_beta(){
    EW_betasq=square(mu_beta)+sigma2_beta;
    //EW_leastSquares =as_scalar(yty-2*ytX*mu_beta +accu(XtX % (Sigma_beta+mu_beta*trans(mu_beta))));
    EW_leastSquares =as_scalar(yty-2*ytX*mu_beta +accu(diagXtX % sigma2_beta) + trans(mu_beta)*XtX*mu_beta);
      // auto time_beta = get_time::now() - start_beta;
      // if(verbose) Rcout<<"Time required:"<<std::chrono::duration_cast<std::chrono::milliseconds>(time_beta).count()<<" ms "<<endl;
  }

  //function to update expected values involving tau
  void update_exp_tau(){
    EW_tau=alpha_tau/beta_tau;
  }

  //function to update expected values involving gamma
  void update_exp_gamma(){
    EW_gamma=alpha_gamma/beta_gamma;
  }


  //function to calculate ELBO
  void calculate_ELBO(){
    // if(verbose) Rcout<<"Calculating ELB.."<<endl;
    // auto start_ELB=get_time::now();

    double ELB_old = ELB;

    vec lgamma_alpha_gamma(g);
    vec digamma_alpha_gamma(g);
    for(int i =0; i<g;i++){
      lgamma_alpha_gamma(i)=lgamma((alpha_gamma(i)));                     // important to directly use log gamma to avoid numerical overflow
      digamma_alpha_gamma(i)=boost::math::digamma(alpha_gamma(i));
    }

    //expected values required in addition (log of Gamma r.v.)
    double EW_logtau = boost::math::digamma(alpha_tau)-log(beta_tau);
    vec EW_loggamma = digamma_alpha_gamma-log(beta_gamma);

    //to get EW_loggamma[annot] and EW_gamma[annot]
    vec EW_loggamma_annot(p);
    vec EW_gamma_annot(p);
    for(int i = 0; i< p; i++){
      int k = annot(i)-1;                          // minus one as annot stars counting at 1 instead of 0
      EW_loggamma_annot(i) = EW_loggamma(k);
      EW_gamma_annot(i) = EW_gamma(k);
    }

    //expectation under variational density of log joint distribution
    double exp_logcondDy=n/2*EW_logtau -0.5*EW_tau*EW_leastSquares-n/2*log(2*M_PI);
    double exp_logcondDbeta=accu(0.5*EW_loggamma_annot-0.5*EW_gamma_annot%EW_betasq)-p/2*log(2*M_PI);
    double exp_logDgamma=accu((r_gamma_mult-1)%EW_loggamma-d_gamma * EW_gamma)-accu(lgamma(r_gamma_mult))+accu(r_gamma_mult*log(d_gamma));
    double exp_logDtau=(r_tau-1)*EW_logtau-d_tau* EW_tau-lgamma(r_tau)+r_tau*log(d_tau);
    double exp_Djoint=exp_logcondDy+exp_logcondDbeta+exp_logDgamma+exp_logDtau;


    //entropy of variational distribution
    double logdet_Sigma = accu(log(sigma2_beta));     //need to calculate sum of log instead of log of product to avoid numerical overflow
    double entropy_beta=p/2*(log(2*M_PI)+1)+0.5*logdet_Sigma;
    double entropy_gamma=accu(alpha_gamma-log(beta_gamma)+ (lgamma_alpha_gamma)+(1-alpha_gamma)%digamma_alpha_gamma); //replace log(tgamma) by lgamma to avoid numeric issues of Inf
    double entropy_tau=alpha_tau-log(beta_tau)+lgamma(alpha_tau)+(1-alpha_tau)*boost::math::digamma(alpha_tau);

    //evidence lower bound
    ELB=exp_Djoint+entropy_beta +entropy_gamma+entropy_tau;
    diff=ELB-ELB_old;

    if(verbose){
      Rcout<<"ELB="<<ELB<<" diff="<<diff<<endl;
    }

  }

};
