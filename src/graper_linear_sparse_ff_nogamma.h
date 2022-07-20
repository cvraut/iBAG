#include <armadillo>
#include  <boost/math/special_functions/digamma.hpp>
#include  <boost/math/special_functions/beta.hpp>

# include <iostream>
# include <chrono>
using get_time = std::chrono::steady_clock ;

using  namespace  arma;
using namespace Rcpp;


class  graper_sparse_ff_nogamma {

private:
  // values remaining constant
  mat  X, XtX;
  vec y, Xty, diagXtX;
  rowvec ytX;
  Row<int> annot;
  double yty;
  int p,n,g;
  vec NoPerGroup;
  double d_tau, r_tau, d_gamma, r_gamma, d_pi, r_pi;
  int max_iter;
  double th;
  bool calcELB, verbose;
  int freqELB;

  // changing values
  vec mu_tildebeta_1, sigma2_tildebeta_0, sigma2_tildebeta_1;         //parameters of the two mixture components in the posterior, mu_tildebeta_0 is zero
  vec psi;                                                                         // probability of s=1
  double EW_tau, ELB;
  double alpha_tau, beta_tau;
  double  alpha_gamma, beta_gamma;
  vec  alpha_pi, beta_pi;
  vec EW_pi;
  vec sigma_beta;
  vec mu_beta;                                                                      // expected values of beta=tildebeta*s
  vec mu_betasq;                                                                   // expected values of beta^2=tildebeta^2*s
  mat Sigma_beta;                                                                   // variance of beta=tildebeta*s
  double EW_gamma;
  vec EW_betatilde;
  double diff;
  int n_iter;
  vec EW_betatildesq;
  double EW_leastSquares;
  vec EW_logfrac_pi;
  vec ELB_trace;
public:

  //initaliser list
  graper_sparse_ff_nogamma(mat X, vec y, Row<int> annot, int g, vec NoPerGroup, double d_tau, double r_tau,
       double d_gamma, double r_gamma, double r_pi, double d_pi, int max_iter, double th, bool calcELB,
                  bool verbose, int freqELB, vec mu_init, vec psi_init):
  X(X)                                  // design matrix
    , XtX(trans(X)*X)
  , y(y)                                // response vector
    , Xty(trans(X)*y)
    , diagXtX(XtX.diag())
    , ytX(trans(y)*X)
    , annot(annot)                        // assignement of each feautre to a group
    , yty(as_scalar(trans(y)*y))
    , p(X.n_cols)                         //number of samples
    , n(X.n_rows)                         //number of samples
    , g(g)                                 // number of groups
    , NoPerGroup(NoPerGroup)               //number of features per group
  , d_tau(d_tau)                        // hyperparameters of gamma distribution for tau
  , r_tau(r_tau)                        // hyperparameters of gamma distribution for tau
  , d_gamma(d_gamma)                    // hyperparameters of gamma distribution for gamma
  , r_gamma(r_gamma*p)                    // hyperparameters of gamma distribution for gamma
  , d_pi(d_pi)                          // hyperparameters of Beta distribution for pi
  , r_pi(r_pi)                          // hyperparameters of Beta distribution for pi
  , max_iter(max_iter)                  // maximal number of iterations
  , th(th)                              //threshold for ELBO to stop iterations
  , calcELB(calcELB)                    //whether to calculate ELBO
  , verbose(verbose)                    //whether to print intermediate messages
  , freqELB(freqELB)                    // freuqency of ELB calculation: each freqELB-th iteration ELBO is calculated
    , mu_tildebeta_1(p)
    , sigma2_tildebeta_0(p)
    , sigma2_tildebeta_1(p)
    , psi(psi_init)                              // is idnetcal to EW_S as Bernoulli
  , EW_tau(r_tau/d_tau)                 //initialise by expected value of a gamma distribution
, ELB(-std::numeric_limits<double>::infinity())                           //evidence lower bound
    , alpha_tau(r_tau+n/2)                //parameter of gamma distribution for tau (stays same in each iteration)
    , alpha_gamma(r_gamma+p/2)             //parameter of gamma distribution for gamma (stays same in each iteration)
    , beta_gamma(d_gamma)
    , alpha_pi(g)
    , beta_pi(g)
    , EW_pi(g)
  , mu_beta(mu_init)                    //initialised randomly
  , EW_gamma(r_gamma/d_gamma)            //initialise by expected value of a gamma distribution, one value per group
  , diff(th+1)                          // to ensure it is larger than th at the beginning
  , n_iter(0)                           // counter of iterations
  , ELB_trace(max_iter)
  {

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

    List results=List::create(Named("EW_beta")=mu_beta,Named("EW_s")=psi, Named("EW_gamma")=EW_gamma,Named("EW_tau")=EW_tau, Named("ELB")=ELB,
                                    Named("alpha_gamma")=alpha_gamma, Named("alpha_tau")=alpha_tau, Named("beta_tau")=beta_tau,
                                    Named("beta_gamma")=beta_gamma, Named("Sigma_beta")=Sigma_beta, Named("EW_pi")=EW_pi,
                                Named("ELB_trace")=ELB_trace);

    return(results);
  }

  //function to do one iteration
  void iterate(){
    n_iter=n_iter+1;                          //increasing counter by 1
    if(verbose) Rcout << "iteration " << n_iter << endl;

    update_param_pi();
    update_exp_pi();
    update_param_beta();
    update_exp_beta();
    update_param_tau();
    update_exp_tau();
    update_param_gamma();
    update_exp_gamma();

    //optional: calculate ELB every freqELB-th step
    if(calcELB && (n_iter%freqELB==0)) calculate_ELBO();
    ELB_trace(n_iter-1)=ELB;
  }


  //function to calculate updated parameters for beta variational distirbution
  void update_param_beta(){
    if(verbose) Rcout << "Updating beta.." << endl;
    auto start_beta=get_time::now();

    vec gamma_annot(p);
    for(int i = 0; i< p; i++) {
      gamma_annot(i)=EW_gamma;      // minus one as annot starts counting at 1 instead of 0
    }
      
    vec EW_logfrac_pi_annot(p);
    for(int i = 0; i< p; i++){
      int k = annot[i]-1;           //group of feautre i
      EW_logfrac_pi_annot(i)=EW_logfrac_pi[k];
    }

    //parameter of normal distribution given s=0
    sigma2_tildebeta_0=1/gamma_annot;
    //parameter of normal distribution given s=1 and probability of s=1
    sigma2_tildebeta_1=1/(EW_tau*diagXtX+gamma_annot);
    
    vec vec1 = X*mu_beta;
    for(int i = 0; i< p; i++){
      //mean of component for s=1
        //  mu_tildebeta_1(i)= sigma2_tildebeta_1(i)* EW_tau * (Xty(i)- accu(XtX.row(i)%trans(mu_beta)) + XtX(i,i)*mu_beta(i));
        // keep track of old mu for efficient update of vec1
        double old_mu_i = mu_beta(i);
        
        mu_tildebeta_1(i)= sigma2_tildebeta_1(i)* EW_tau * (Xty(i)- accu(X.col(i) % vec1) + diagXtX(i)*mu_beta(i));

        
        //probability of s=1
        double term= EW_logfrac_pi_annot(i) + 0.5*log(sigma2_tildebeta_1(i)) - 0.5*log(sigma2_tildebeta_0(i)) + 0.5*mu_tildebeta_1(i)*mu_tildebeta_1(i)*(1/sigma2_tildebeta_1(i));
        psi(i) = 1/(1+exp(-term));
        //expected values of beta=tildebeta*s
        mu_beta(i) = mu_tildebeta_1(i)*psi(i);
        
    //update vec1 (only in new coordinate of mu to avoid recomputing the full pxp product and get linear complexity)
    vec1 = vec1 + (mu_beta(i) - old_mu_i)*X.col(i);

    }

    auto time_beta = get_time::now() - start_beta;
    if(verbose) Rcout<<"Time required:"<<std::chrono::duration_cast<std::chrono::milliseconds>(time_beta).count()<<" ms "<<endl;
  }

  //function to calculate updated parameters for tau variational distribution
  void update_param_tau(){
    beta_tau=d_tau+0.5*EW_leastSquares;
  }

  //function to calculate updated parameters for gamma variational distribution
  void update_param_gamma(){
    beta_gamma = d_gamma;

    for(int i = 0; i< p; i++){
      beta_gamma=beta_gamma+0.5*EW_betatildesq[i];
    }
  }

  //function to calculate updated parameters for pi variational distribution
  void update_param_pi(){
    alpha_pi.fill(d_pi);
    beta_pi.fill(r_pi);

    for(int i = 0; i< p; i++){
      int k = annot[i]-1;                          // minus one as annot stars counting at 1 instead of 0
      alpha_pi[k]=alpha_pi[k]+psi[i];
      beta_pi[k]=beta_pi[k]+(1-psi[i]);
    }
  }

  //function to update expected values of beta
  void update_exp_beta(){
      if(verbose) Rcout << "Updating expected values containing beta.." << endl;
      auto start_beta=get_time::now();
    EW_betatildesq=psi%(square(mu_tildebeta_1) +sigma2_tildebeta_1) + (1-psi)% (sigma2_tildebeta_0);


    //expected value of beta^2= tildebeta^2*s
    mu_betasq=psi%(square(mu_tildebeta_1)+sigma2_tildebeta_1);

    //variance of beta=tildebeta*s
    Sigma_beta = speye(p,p);
    Sigma_beta.diag()=(mu_betasq-square(mu_beta));

    //expected value of least squares expression
      //EW_leastSquares =as_scalar(yty-2*ytX*mu_beta +accu(XtX % (Sigma_beta+mu_beta*trans(mu_beta))));
      EW_leastSquares =as_scalar(yty-2*ytX*mu_beta +accu(XtX % Sigma_beta) + trans(mu_beta)*XtX*mu_beta);
      auto time_beta = get_time::now() - start_beta;
      if(verbose) Rcout<<"Time required:"<<std::chrono::duration_cast<std::chrono::milliseconds>(time_beta).count()<<" ms "<<endl;
  }

  //function to update expected values of tau
  void update_exp_tau(){
    EW_tau=alpha_tau/beta_tau;
  }

  //function to update expected values of gamma
  void update_exp_gamma(){
    EW_gamma=alpha_gamma/beta_gamma;
  }

  //function to update expected values of gamma
  void update_exp_pi(){
    vec digamma_alpha_pi(g);
    vec digamma_beta_pi(g);
    for(int i =0; i<g;i++){
      digamma_alpha_pi(i)=boost::math::digamma(alpha_pi(i));                     // important to directly use log gamma to avoid numerical overflow
      digamma_beta_pi(i)=boost::math::digamma(beta_pi(i));
    }
    EW_logfrac_pi=digamma_alpha_pi-digamma_beta_pi;                             //\mathbb{E}\log\frac{\pi_k}{1-\pi_k} 
    EW_pi = alpha_pi/(alpha_pi+beta_pi);
  }


  //function to calculate ELBO
  void calculate_ELBO(){
     if(verbose) Rcout<<"Calculating ELB.."<<endl;
     auto start_ELB=get_time::now();

    double ELB_old = ELB;

    double lgamma_alpha_gamma;
    double digamma_alpha_gamma;
    digamma_alpha_gamma = boost::math::digamma(alpha_gamma);
    lgamma_alpha_gamma=lgamma((alpha_gamma));                     // important to directly use log gamma to avoid numerical overflow

      
    vec digamma_alpha_pi(g);
    vec digamma_beta_pi(g);
    vec digamma_alphaplusbeta_pi(g);
    vec logbeta_alpha_beta_pi(g);
    for(int i =0; i<g;i++){
      digamma_alpha_pi(i)=boost::math::digamma(alpha_pi(i));
      digamma_beta_pi(i)=boost::math::digamma(beta_pi(i));
      digamma_alphaplusbeta_pi(i)=boost::math::digamma(alpha_pi(i)+beta_pi(i));
      //logbeta_alpha_beta_pi(i)=log(boost::math::beta(alpha_pi(i),beta_pi(i)));
      logbeta_alpha_beta_pi(i)=lgamma(alpha_pi(i)) + lgamma(beta_pi(i)) - lgamma(alpha_pi(i)+beta_pi(i));        //try whether there is also numerical overflow

     }

    //expected values required in addition (log of Gamma r.v.)
    //Note: beta_gamma and beta_tau are same for all groups and constant across iterations
    double EW_logtau = boost::math::digamma(alpha_tau)-log(beta_tau);
    double EW_loggamma = digamma_alpha_gamma-log(beta_gamma);
    vec EW_logpi = digamma_alpha_pi-digamma_alphaplusbeta_pi;
    vec EW_logOneMinusPi= digamma_beta_pi-digamma_alphaplusbeta_pi;

    //to get EW_loggamma[annot] and EW_gamma[annot]
    vec EW_loggamma_annot(p);
    vec EW_gamma_annot(p);
    vec EW_logpi_annot(p);
    vec EW_logOneMinusPi_annot(p);
    for(int i = 0; i< p; i++){
      int k = annot(i)-1;                          // minus one as annot stars counting at 1 instead of 0
      EW_loggamma_annot(i) = EW_loggamma;
      EW_gamma_annot(i) = EW_gamma;
      EW_logpi_annot(i) = EW_logpi(k);
      EW_logOneMinusPi_annot(i) = EW_logOneMinusPi(k);

    }

    //expectation under variational density of log joint distribution
    double exp_logcondDy = n/2*EW_logtau -0.5*EW_tau*EW_leastSquares-n/2*log(2*M_PI);
    double exp_logcondDbeta = accu(0.5*EW_loggamma_annot-0.5*EW_gamma_annot%EW_betatildesq)-p/2*log(2*M_PI);
    double exp_logcondDs = accu(psi%EW_logpi_annot+(1-psi)%EW_logOneMinusPi_annot);
    double exp_logDgamma = accu((r_gamma-1)*EW_loggamma-d_gamma * EW_gamma)-g*lgamma(r_gamma)+g*r_gamma*log(d_gamma);
    double exp_logDtau = (r_tau-1)*EW_logtau-d_tau* EW_tau-lgamma(r_tau)+r_tau*log(d_tau);
    //double exp_logDpi = accu((d_pi-1)*EW_logpi +(r_pi-1)*EW_logOneMinusPi)-g*log(boost::math::beta(d_pi,r_pi));
    double exp_logDpi = accu((d_pi-1)*EW_logpi +(r_pi-1)*EW_logOneMinusPi)-g*(lgamma(d_pi)+lgamma(r_pi)-lgamma(d_pi+r_pi));
    double exp_Djoint = exp_logcondDy+exp_logcondDbeta+exp_logcondDs+exp_logDgamma+exp_logDtau+exp_logDpi;

    //entropy of variational distribution
    vec entropy_s_indiv = psi%log(psi) +(1-psi)%log(1-psi);
    //entropy of bernoulli s becomes nan for psi=0 or 1, here the limit it zero in both cases
    for(int i = 0; i< p; i++){
      if((psi(i)==1) || (psi(i)==0)) entropy_s_indiv(i) = 0;
    }

    double entropy_tildebeta_s = p/2*(log(2*M_PI)+1)+0.5*accu(psi%log(sigma2_tildebeta_1) +(1-psi)%log(sigma2_tildebeta_0))-accu(entropy_s_indiv);
    double entropy_gamma = accu(alpha_gamma-log(beta_gamma)+(lgamma_alpha_gamma)+(1-alpha_gamma)*digamma_alpha_gamma); //replace log(tgamma) by lgamma to avoid numeric issues of Inf
    double entropy_tau = alpha_tau-log(beta_tau)+lgamma(alpha_tau)+(1-alpha_tau)*boost::math::digamma(alpha_tau);
    double entropy_pi = accu(logbeta_alpha_beta_pi - (alpha_pi-1)%digamma_alpha_pi -
                             (beta_pi-1)%digamma_beta_pi +(alpha_pi+beta_pi-2)%digamma_alphaplusbeta_pi);

    //evidence lower bound
    ELB=exp_Djoint + entropy_tildebeta_s + entropy_gamma + entropy_tau + entropy_pi;
    diff=ELB-ELB_old;
    //
    auto time_ELB = get_time::now() - start_ELB;
    if(verbose) Rcout<<"Time required:"<<std::chrono::duration_cast<std::chrono::milliseconds>(time_ELB).count()<<" ms "<<endl;

    if(verbose){
      Rcout<<"ELB="<<ELB<<endl;
      Rcout<<"ELB improved by "<<diff<<endl;
      Rcout<<endl;
    }

  }

};
