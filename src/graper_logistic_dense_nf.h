#include <armadillo>
#include  <boost/math/special_functions/digamma.hpp>

# include <iostream>
# include <chrono>
using get_time = std::chrono::steady_clock ;

using  namespace  arma;
using namespace Rcpp;


class  graper_logistic_nf {

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
  bool calcELB, verbose;
  int freqELB;
  List ListOfOuterX;
  vec term4betamu;

  // changing values
  double ELB;
  vec  alpha_gamma, beta_gamma;
  vec xi;
  mat Sigma_beta;
  vec mu_beta;
  vec EW_gamma;
  double diff;
  int n_iter;
  vec EW_betasq;
  vec ELB_trace;
public:

  //initaliser list
  graper_logistic_nf(mat X, vec y, Row<int> annot, int g, vec NoPerGroup,
       double d_gamma =0.001, double r_gamma =0.001, int max_iter=1000, double th=1e-7, bool calcELB=true, bool verbose=true,
       int freqELB =10):
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
  , freqELB(freqELB)                    // freuqency of ELB calculation: each freqELB-th iteration ELBO is calculated
  , ListOfOuterX(n)
  , term4betamu(p)
  , ELB(-std::numeric_limits<double>::infinity())                           //evidence lower bound
    , alpha_gamma(g)                      //parameter of gamma distribution for tau (stays same in each iteration)
  , beta_gamma(g)
    , xi(n)                               // variational parameter, initialised to 0, should better be close to yX\beta
  , EW_gamma(g)                         //initialise by expected value of a gamma distribution, one value per group
  , diff(th+1)                          // to ensure it is larger than th at the beginning
  , n_iter(0)                           // counter of iterations
  , ELB_trace(max_iter)
  {
    EW_gamma.fill(r_gamma/d_gamma);
    alpha_gamma=r_gamma+NoPerGroup/2;
    term4betamu.fill(0);
    xi.fill(0);
    //calculate often used quantities - slow
    for(int i = 0; i< n; i++) {
      ListOfOuterX(i) = X.row(i).t()*X.row(i);
      term4betamu = term4betamu + (y(i)-0.5)*X.row(i).t();
      }
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

    List results=List::create(Named("EW_beta")=mu_beta, Named("EW_gamma")=EW_gamma, Named("ELB")=ELB,
                                    Named("alpha_gamma")=alpha_gamma,
                                    Named("beta_gamma")=beta_gamma, Named("Sigma_beta")=Sigma_beta, Named("ELB_trace")=ELB_trace);


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
    update_param_xi();


    //optional: calculate ELB every freqELB-th step
    if(calcELB & (n_iter%freqELB==0)) calculate_ELBO();
    ELB_trace(n_iter-1)=ELB;

  }


  //function to calculate updated parameters for beta variational distirbution
  void update_param_beta(){
    if(verbose) Rcout << "Updating beta.." << endl;
    auto start_beta=get_time::now();

    sp_mat A= speye(p,p);               //create diagonal matrix with entries of EW_gamma according to group membership - essential to use speye INSTEAD OF EYE for speed
    sp_mat Id_n = speye(n,n);
    vec gamma_annot(p);
    for(int i = 0; i< p; i++) {
      gamma_annot(i)=EW_gamma(annot(i)-1);      // minus one as annot starts counting at 1 instead of 0
    }
    mat term1 =zeros(p,p);
    for(int i = 0; i< n; i++) {
      mat outerXi = ListOfOuterX(i);
      term1=term1+lambda(xi(i)) * outerXi;
    }

      A.diag() = gamma_annot;
      Sigma_beta = inv_sympd(2*term1 +A);                  //inverse of symmetric, positive definite matrix -slow


    mu_beta =Sigma_beta*term4betamu;

    auto time_beta = get_time::now() - start_beta;
    if(verbose) Rcout<<"Time required:"<<std::chrono::duration_cast<std::chrono::milliseconds>(time_beta).count()<<" ms "<<endl;
  }

  //function to calculate updated parameters for tau variational distribution
  void update_param_xi(){
    for(int i = 0; i< n; i++) {
      xi(i) = sqrt(as_scalar(X.row(i)*(Sigma_beta + mu_beta*mu_beta.t())*X.row(i).t()));    //slow
    }
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
    EW_betasq=square(mu_beta)+Sigma_beta.diag();
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
