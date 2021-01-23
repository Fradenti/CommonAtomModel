#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// Auxiliary functions -------------------------------------------------------------------

// [[Rcpp::export]]
arma::cube heat_mat(arma::mat x) {
  int J = size(x)[0], NSIM = size(x)[1];
  arma::cube CU(J, J, NSIM);
  
  for(int i=0; i<NSIM ;i++){
    for(int l=0; l<J ;l++){
      for(int k=0; k<J ;k++){
        if(x(l,i)==x(k,i)){CU(l, k, i)=1;
        }else{
          CU(l, k, i)=0; } 
      }}}
  return CU;
}

// [[Rcpp::export]]
arma::vec table_factor_SLICE(arma::colvec Ci, int L){
  arma::vec tab(L);
  for(int l=0; l < L; l++){
    tab[l] =  accu( Ci==(l+1) );
  }
  return(tab);
}

// Discrete Log-likelihood with Library Size -------------------------------------------------------------------
// [[Rcpp::export]]
double log_Likelihood(double yobs, arma::rowvec theta, double gamma){
  
  double res=0.0;
  double current_mu = (theta[0])*gamma;
  double current_sigma = sqrt(theta[1])*gamma, a1,a2;
  

    bool check1 = (yobs+1 - current_mu)/(current_sigma) > 5;
    bool check2 = (yobs   - current_mu)/(current_sigma) > 5;
    
    if((check1+check2)==0){
      
      if(yobs==0){
        res = R::pnorm(1, current_mu, current_sigma, 1, 0);
      }else{
        a1  = R::pnorm(yobs+1, current_mu, current_sigma, 1, 0);
        a2  = R::pnorm(yobs  , current_mu, current_sigma, 1, 0);
        res = a1-a2;}
    }else{
      if(yobs==0){
        res = R::pnorm(1, current_mu, current_sigma, 0, 0);
      }else{
        a1  = R::pnorm(yobs+1, current_mu, current_sigma, 0, 0);
        a2  = R::pnorm(yobs  , current_mu, current_sigma, 0, 0);
        res = a2-a1;
      }}
  return log(res);
}



// Slice sampler functions --------------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec g_cpp(arma::colvec j, double kappa){
  arma::colvec jj = log(1-kappa) + (j-1)*log(kappa);
  return( exp(jj));
}


// Full conditionals --------------------------------------------------------------------------

// Update Zj
// [[Rcpp::export]]
arma::colvec UPD_Zj(arma::colvec Uj, 
                    arma::colvec Uij,
                    arma::colvec xi_z, 
                    arma::colvec xi_c,
                    arma::colvec pi_z,
                    arma::colvec cij,
                    arma::mat omega,
                    arma::colvec y_group,
                    int NN_z, int J){
  arma::colvec newZj(J), p_zj_k(NN_z),
  possible_label = arma::linspace<arma::vec>(1, NN_z, NN_z);
  
  for(int q=0; q<J; q++){
    arma::uvec      ind = find(y_group==(q+1));
    arma::colvec subUij = Uij.elem(ind), 
      subCij = cij.elem(ind);
    arma::uvec indCij   = arma::conv_to<arma::uvec>::from(subCij-1);
    arma::mat subOmega  = omega.rows(indCij);
    
    for(int k=0; k<NN_z; k++){
      p_zj_k[k] =   log( Uj[q] < xi_z[k] ) - log(xi_z[k]) + log(pi_z[k]) +
        accu( log( subUij < xi_c.elem(indCij))) +
        accu(log( subOmega.col(k))) - 
        accu(log(xi_c.elem(indCij)));
    }
    arma::colvec pp = exp(p_zj_k-max(p_zj_k));
    newZj[q] = RcppArmadillo::sample(possible_label, 1, 1, pp)[0];      
  }
  return newZj;
}


// Update Cij
// [[Rcpp::export]]
arma::mat UPD_tsl0_for_cij( arma::colvec y_obser,
                            arma::colvec Uij,
                            arma::colvec xi_c,
                            arma::mat omega,
                            arma::colvec zj_pg,
                            arma::mat tsl0,
                            arma::colvec gamma,
                            int N, int NN_c){
  arma::colvec possible_label = arma::linspace<arma::vec>(1, NN_c, NN_c);
  arma::mat tsl0_tmp(N,2);
  arma::colvec p(NN_c);
  for(int i=0; i<N; i++){
    for(int k=0; k<NN_c; k++){
      p[k] = log(xi_c[k] > Uij[i]) + log(omega(k,zj_pg[i]-1)) +  
        log_Likelihood( y_obser[i],  tsl0.row(k), gamma[i]) - log(xi_c[k]);    
    }
    arma::colvec pp = exp(p-max(p));
    int IND = RcppArmadillo::sample(possible_label, 1, 1, pp)[0];
    tsl0_tmp.row(i) = tsl0.row(IND-1); 
  }
  return(tsl0_tmp);
}




double rt_cpp(double nu, double lambda){
  double TAA;
  TAA = R::rnorm(lambda,1.0) / exp(.5*log( R::rchisq(nu)/nu ));
  return(TAA);
}

// [[Rcpp::export]]
arma::colvec SB_given_u2(arma::colvec V) {
  
  int N = V.size();
  arma::colvec pi(N), mV2(N);
  arma::colvec mV = 1-V;
  mV2=arma::shift(cumprod(mV),+1);
  mV2(0) =1;
  return(V%mV2);
}


// Update tsl0
// [[Rcpp::export]]
arma::mat UPD_tsl0(arma::colvec y_LAT,
                   arma::colvec cij,
                   double a0, double b0, double k0, double m0,
                   int NN_c, int J){
  arma::colvec cij_cpp = cij -1 ;
  arma::mat tsl0(NN_c,2);
  double ybar_i, ss_i;
  double astar, bstar, mustar, kstar;
  for(int i = 0; i<NN_c; i++ ){
    arma::uvec   ind = find(cij_cpp==i);
    arma::colvec YYY = y_LAT.elem(ind);
    int          n_i = YYY.n_elem;
    if(n_i > 0) {  ybar_i = mean(YYY);} else {ybar_i = 0;}
    if(n_i > 1) {  ss_i   = accu( pow((YYY-ybar_i),2) );} else {ss_i = 0;}
    astar = (a0 + n_i / 2);
    bstar = b0 + .5 * ss_i + ((k0 * n_i) * (ybar_i - m0)* (ybar_i - m0))  / ( 2 * (k0+n_i) );
    tsl0(i,1) = 1 / rgamma(1, astar, 1/bstar)[0];
    mustar    = (k0*m0+ybar_i*n_i)/(k0+n_i);
    kstar     = k0 + n_i;
    tsl0(i,0) = rt_cpp(2*astar, 0.0) * sqrt(bstar/(kstar*astar))+mustar;  
  }
  return(tsl0);
}




// Update omega
// [[Rcpp::export]]
arma::mat UPD_omega(arma::colvec cij, 
                    arma::colvec zj_pg,
                    int NN_c, int NN_z, double beta){
  arma::colvec zj_pg_cpp = zj_pg -1 ,
    cij_cpp = cij -1;
  arma::mat v_omega(NN_c, NN_z), omega(NN_c, NN_z);
  v_omega.fill(0); v_omega.fill(0);
  
  for(int jj=0; jj<NN_z; jj++){
    for(int ii=0; ii<NN_c;  ii++){
      v_omega(ii,jj) = R::rbeta( 1 + accu(zj_pg_cpp == jj && cij_cpp == ii),
              beta + accu(zj_pg_cpp == jj && cij_cpp > ii));
    }
    omega.col(jj) = SB_given_u2(v_omega.col(jj));
  }
  return(omega); 
}

// Update pi
// [[Rcpp::export]]
arma::colvec UPD_Pi_v_z(arma::colvec zj, 
                 int NN_z, double alpha){
  arma::colvec v_z(NN_z), pi_z(NN_z);
  for(int j=0; j<NN_z; j++){
    v_z[j] = R::rbeta(1 + accu(zj == (j+1)), alpha + accu(zj > (j+1)) );
  }
  return v_z;
}

// Update latent variable ---------------------------------------------
// [[Rcpp::export]]
arma::colvec Upd_LAT(arma::colvec y_obser,  
                     arma::colvec cij, arma::colvec gamma,
                     arma::mat tsl0, int N){
  
  
  arma::colvec ylat(N), index(N), cij_cpp = cij-1;
  double    current_sigma, current_mu;
  arma::vec thresh(2*N), u(N);
  arma::colvec Tsl01_selected = tsl0.col(0),
               Tsl02_selected = tsl0.col(1);
  
  for(int i=0; i<N; i++){
    current_mu    =      tsl0(cij_cpp[i],0) *gamma[i];
    current_sigma = sqrt(tsl0(cij_cpp[i],1))*gamma[i]; 
    
    bool check1 = (y_obser[i]+1 - current_mu)/(current_sigma) > 5;
    bool check2 = (y_obser[i]   - current_mu)/(current_sigma) > 5;
    
    if((check1+check2)==0){
      thresh[N+i]= R::pnorm(y_obser[i]+1, current_mu, current_sigma, 1, 1);
      thresh[i]  = R::pnorm(y_obser[i]  , current_mu, current_sigma, 1, 1);
      u[i]       = R::runif(thresh[i], thresh[N+i]);
      ylat[i]    = R::qnorm(u[i], current_mu, current_sigma,1,1);
    }else{
      thresh[N+i]= R::pnorm(y_obser[i]+1, current_mu, current_sigma, 0, 1);
      thresh[i]  = R::pnorm(y_obser[i]  , current_mu, current_sigma, 0, 1);
      u[i]       = R::runif(thresh[N+i], thresh[i]);
      ylat[i]    = R::qnorm(u[i], current_mu, current_sigma,0,1);
    }
  }
  return ylat;
}



// Update labels - collapsed Gibbs

// [[Rcpp::export]]
arma::colvec UPD_Zj_collapsed_uij(arma::colvec Uj, 
                                  arma::colvec xi_z, 
                                  arma::colvec xi_c,
                                  arma::colvec pi_z,
                                  arma::colvec cij,
                                  arma::mat omega,
                                  arma::colvec y_group,
                                  int NN_z, int J){
  arma::colvec newZj(J), p_zj_k(NN_z),
  possible_label = arma::linspace<arma::vec>(1, NN_z, NN_z);
  
  for(int q=0; q<J; q++){
    arma::uvec      ind = find(y_group==(q+1));
    arma::colvec    subCij = cij.elem(ind);
    arma::uvec indCij   = arma::conv_to<arma::uvec>::from(subCij-1);
    arma::mat subOmega  = omega.rows(indCij);
    
    for(int k=0; k<NN_z; k++){
      p_zj_k[k] =   log( Uj[q] < xi_z[k] ) - log(xi_z[k]) + log(pi_z[k]) +
        accu(log( subOmega.col(k))) ;
    }
    arma::colvec pp = exp(p_zj_k-max(p_zj_k));
    newZj[q] = RcppArmadillo::sample(possible_label, 1, 1, pp)[0];      
  }
  return newZj;
}

