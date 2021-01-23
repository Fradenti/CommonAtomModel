int rintnunif_log(arma::vec lweights, int NN_z){
  
  double u = arma::randu();
  arma::vec probs(NN_z);
  
  for(int k = 0; k < NN_z; k++) {
    probs(k) = 1 / sum(exp(lweights - lweights(k)));
  }
  
  for(int k = 1; k <= NN_z; k++) {
    if(u <= probs[k-1]) {
      return k;
    }
  }
}

// Update g
// [[Rcpp::export]]
arma::colvec g_cpp(arma::colvec j, double kappa){
  arma::colvec jj = log(1-kappa) + (j-1)*log(kappa);
  return( exp(jj) );
}



// [[Rcpp::export]]
arma::mat per_Update_omega_matpesi(arma::colvec cij, 
                                   arma::colvec zj_pg,
                                   int NN_c, int NN_z, double beta){
  arma::colvec zj_pg_cpp = zj_pg -1 ,
    cij_cpp = cij -1;
  arma::mat m_cij(NN_c, NN_z);
  m_cij.fill(0);
  
  for(int jj=0; jj<NN_z; jj++){
    for(int ii=0; ii<NN_c;  ii++){
      m_cij(ii,jj) = accu(zj_pg_cpp == jj && cij_cpp == ii);
    }
    
  }
  return(m_cij); 
}


// Update pi
// [[Rcpp::export]]
arma::mat Update_Pi(arma::colvec zj, 
                    int NN_z, double alpha){
  arma::colvec v_z(NN_z), pi_z(NN_z);
  for(int j=0; j<NN_z; j++){
    v_z[j] = R::rbeta(1 + accu(zj == (j+1)), alpha + accu(zj > (j+1)) );
  }
  pi_z= SB_given_u2(v_z);
  return pi_z;
}




/*
 
 // [[Rcpp::export]]
 double Update_alpha_precision_EscWest(double N, double K, double a, double b, double old_aplha){
 
 
 double logeta = log( R::rbeta(old_alpha + 1,N) );  
 double Q      = (a_+K-1)/( N * (b -logeta) );
 double pi_eta = Q/(1+Q);
 double u = R::runif(0,1);
 if()------
 
 alpha  <- if(u<pi_eta){  
 rgamma(n = 1, shape = a_alpha+k,  rate = b_alpha-logeta)
 }else{ 
 rgamma(n = 1, shape = a_alpha+k-1,rate = b_alpha-logeta) 
 }
 
 
 
 return(alpha);  
 }
 */







// Update Cij
// [[Rcpp::export]]
arma::mat Update_theta_for_cij( arma::colvec y_obser,
                                arma::colvec Uij,
                                arma::colvec xi_c,
                                arma::mat omega,
                                arma::colvec zj_pg,
                                arma::mat theta,
                                int N, int NN_c, arma::colvec possible_label){
  arma::mat theta_tmp(N,2);
  arma::colvec p(NN_c);
  
  for(int i=0; i<N; i++){
    
    for(int k=0; k<NN_c; k++){
      
      p[k] = log(xi_c[k] > Uij[i]) + log(omega(k,zj_pg[i]-1)) +  
        R::dnorm(y_obser[i], theta(k,0), sqrt(theta(k,1)), 1) - log(xi_c[k]);    
      
    }
    
    if(arma::is_finite(max(p))){
      arma::colvec pp = exp(p-max(p));
      int IND = RcppArmadillo::sample(possible_label, 1, 1, pp)[0];
      theta_tmp.row(i) = theta.row(IND-1); 
    }else{
      int IND = RcppArmadillo::sample(possible_label, 1, 1)[0];
      theta_tmp.row(i) = theta.row(IND-1); 
    }
  }
  
  return(theta_tmp);
}