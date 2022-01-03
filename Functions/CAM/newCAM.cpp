#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


///////////////////////////////////////////////////// Useful
// [[Rcpp::export]]
arma::mat PSM(arma::mat inds){
  int nsim = inds.n_rows;
  int n= inds.n_cols;
  arma::mat PSM(n,n), D(n,n); 
  D.eye(); PSM.zeros();
  
  for(int i=0; i<n; i++){
    for(int j=i+1; j<n; j++){
    arma::colvec Z = (inds.col(i)-inds.col(j));
    arma::uvec success = find(Z==0);
    PSM(i,j) = success.n_elem;
    PSM(j,i) = PSM(i,j);
    }
  }
  
  return(PSM/nsim + D);
  
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

/////////////////////////////////////////////////////////// Update pi
// [[Rcpp::export]]
arma::colvec Update_Distributional_Sticks(arma::colvec zj, 
                        int NN_z, double alpha){
  arma::colvec v_z(NN_z);
  for(int j=0; j<(NN_z); j++){
    v_z[j] = R::rbeta(1 + accu(zj == (j+1)), alpha + accu(zj > (j+1)) );
  }
  //v_z[NN_z-1] =1.; memento: this is wrong
  return v_z; 
}

/////////////////////////////////////////////////////////// Update Omega
// [[Rcpp::export]]
arma::mat Update_Observational_Sticks(arma::colvec cij, 
                       arma::colvec zj_pg,
                       int NN_c, 
                       int NN_z, 
                       double beta){
  arma::colvec zj_pg_cpp = zj_pg -1 ,  cij_cpp = cij -1;
  
  arma::mat v_omega(NN_c, NN_z);
  
  v_omega.fill(0); v_omega.fill(0);
  
  for(int jj=0; jj<NN_z; jj++){
    for(int ii=0; ii<NN_c;  ii++){
      v_omega(ii,jj) = 
          R::rbeta( 1 + accu(zj_pg_cpp == jj && cij_cpp == ii),
                 beta + accu(zj_pg_cpp == jj && cij_cpp >  ii));
    }
   // v_omega(NN_c-1,jj) = 1.; // new line
   // omega.col(jj)      = SB_given_u2(v_omega.col(jj));
  }
  return(v_omega); 
}



// [[Rcpp::export]]
arma::mat Update_omega(arma::colvec cij, 
                       arma::colvec zj_pg,
                       int NN_c, 
                       int NN_z, 
                       double beta){
  arma::colvec zj_pg_cpp = zj_pg -1 ,  cij_cpp = cij -1;
  
  arma::mat v_omega(NN_c, NN_z), omega(NN_c, NN_z);
  
  v_omega.fill(0); v_omega.fill(0);
  
  for(int jj=0; jj<NN_z; jj++){
    for(int ii=0; ii<NN_c;  ii++){
      v_omega(ii,jj) = 
        R::rbeta( 1 + accu(zj_pg_cpp == jj && cij_cpp == ii),
                  beta + accu(zj_pg_cpp == jj && cij_cpp >  ii));
    }
    // v_omega(NN_c-1,jj) = 1.; // new line
    omega.col(jj)      = SB_given_u2(v_omega.col(jj));
  }
  return(omega); 
}


//////////////////////////////////////////////////////////// Update Zj
// [[Rcpp::export]]
arma::colvec Update_Zj_v2(arma::colvec Uj,    // uij collapsed
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
    newZj[q] =   RcppArmadillo::sample(possible_label, 1, 1, pp)[0];      
  }
  return newZj;
}

/*
 // Update Zj
 // [[Rcpp::export]]
 arma::colvec Update_Zj(arma::colvec Uj, 
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
 */


////////////////////////////////////////////////////////////////////// // Update Cij
// [[Rcpp::export]]
arma::colvec Update_Cij( arma::colvec y_obser,
                               arma::colvec Uij,
                               arma::colvec xi_c,
                               arma::mat omega,
                               arma::colvec zj_pg,
                               arma::mat theta,
                               int N, int NN_c){
  arma::colvec possible_label = arma::linspace<arma::vec>(1, NN_c, NN_c);
  arma::colvec IND(N);
  arma::colvec p(NN_c);
  
  for(int i=0; i<N; i++){
    for(int k=0; k<NN_c; k++){
      p[k] = 
        log(xi_c[k] > Uij[i]) + 
        log(omega(k,zj_pg[i]-1)) +  
        R::dnorm(y_obser[i], theta(k,0), sqrt(theta(k,1)), 1) - 
        log(xi_c[k]);    
    }
    
    if(arma::is_finite(max(p))){
      arma::colvec pp = exp(p-max(p));
      IND(i) = RcppArmadillo::sample(possible_label, 1, 1, pp)[0];
    }else{
      IND(i) = RcppArmadillo::sample(possible_label, 1, 1)[0];
    }
    
  }
  
  return(IND);
}




////////////////////////////////////////////////////////////////////  Update theta
// [[Rcpp::export]]
arma::mat Update_theta(arma::colvec y_obser,
                   arma::colvec cij,
                   double a0, double b0, double k0, double m0,
                   int NN_c, int J){
  arma::colvec cij_cpp = cij -1 ;
  arma::mat theta(NN_c,2);
  
  double ybar_i, ss_i;
  double astar, bstar, mustar, kstar;
  
  for(int i = 0; i<NN_c; i++ ){
    arma::uvec   ind = find(cij_cpp==i);
    arma::colvec YYY = y_obser.elem(ind);
    
    int          n_i = YYY.n_elem;
    if(n_i > 0) {  ybar_i = mean(YYY);} else {ybar_i = 0;}
    
    if(n_i > 1) {  ss_i   = accu( pow((YYY-ybar_i),2) );} else {ss_i = 0;}
    
      astar = (a0 + n_i / 2);
      bstar = b0 + .5 * ss_i + ((k0 * n_i) * (ybar_i - m0)* (ybar_i - m0))  / ( 2 * (k0+n_i) );
      theta(i,1) = 1 / rgamma(1, astar, 1/bstar)[0];
      mustar     = (k0*m0+ybar_i*n_i)/(k0+n_i);
      kstar      = k0 + n_i;
      theta(i,0) = rt_cpp(2*astar, 0.0) * sqrt(bstar/(kstar*astar))+mustar;  
    }
  return(theta);
}





















/*
  // Update Cij
  // [[Rcpp::export]]
  arma::colvec Update_ZjCij( arma::colvec y_obser,
                        arma::colvec y_group,
                        arma::colvec Uj,  
                        arma::colvec Uij,
                        arma::colvec xi_z, 
                        arma::colvec xi_c,
                        arma::mat omega,
                        arma::mat theta,
                        int N, int NN_c,
                        arma::colvec pi_z,
                        int NN_z, int J){
    int all = NN_c * NN_z ;
    arma::colvec possible_label = arma::linspace<arma::vec>(1, all, all);
    arma::colvec IND(N);
    arma::mat p(NN_c,NN_z);
    arma::colvec Mij(N), Zj(J);
    
    
    
    for(int i=0; i<N; i++){
      int q = y_group[i]-1;
      for(int l=0; l<NN_z; l++){
        
        double Prob_L = log(pi_z[l]) + log( xi_z[l] >  Uj[q] ) - log(xi_z[l]) ;
        
        for(int k=0; k<NN_c; k++){
          
          double log_prob_Normal = R::dnorm(y_obser[i], theta(k,0), sqrt(theta(k,1)), 1);
          
          
          p(k,l) = 
            log( xi_c[k] > Uij[i] ) +   log(omega(k,l))  - 
            log(xi_c[k]) + Prob_L + log_prob_Normal;
        }
      }
      arma::colvec longp = p.as_col();
      
      
      if(arma::is_finite(max(longp))){
        arma::colvec pp = exp(longp-max(longp));
        IND(i) = RcppArmadillo::sample(possible_label, 1, 1, pp)[0];
      }else{
        IND(i) = RcppArmadillo::sample(possible_label, 1, 1)[0];
      }
      
    }
    
    return(IND);
  }
*/

