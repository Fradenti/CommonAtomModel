#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec Likelihood_row(NumericVector yobs, arma::rowvec theta){
 return  dnorm4(yobs,theta[0],sqrt(theta[1]),0);
}

// [[Rcpp::export]]
arma::colvec Likelihood_col(NumericVector yobs, arma::vec theta){
  return  dnorm4(yobs,theta[0],sqrt(theta[1]),0);
}

// [[Rcpp::export]]
arma::rowvec log_Likelihood(NumericVector yobs, arma::vec theta){
  return  dnorm4(yobs,theta[0],sqrt(theta[1]),1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec UPD_DISTR_LAB(arma::colvec y_obser, 
                        arma::colvec y_group,
                        arma::cube tslk, 
                        arma::mat tsl0, 
                        arma::vec psk, 
                        arma::colvec rho,   
                        arma::mat oslk,  
                        arma::vec osl0, 
                        int J, 
                        int K, 
                        int L) {
  
  arma::vec actual_label(J), 
            possible_label = arma::linspace<arma::vec>(1, K, K),
            prob1(K), prob2(K), prob3(K), prob4(K);

  for(int j=0; j<J; j++){ // for each individual (group), we need to decide in which cluster he is allocated
    
    arma::uvec inds      = find( y_group == (j+1) ) ;
    arma::vec  yj        = y_obser.elem(inds); // take all Y in group j
    int n_yj             = size(yj)[0];
    
    // container
    arma::mat Lik_forall_L(L,n_yj), 
              Lik_forall_L0(L,n_yj) , 
              OSL(L, n_yj),       
              OSL0(L, n_yj), 
              WtimesLik(L, n_yj), 
              WtimesLik0(L, n_yj);
    
    OSL0 = arma::repelem(osl0, 1, n_yj);
    
    for(int k=0; k<K; k++){
      
      OSL = arma::repelem(oslk.col(k), 1, n_yj); // matrix with SB weights (one for each atom till L), repeated by column, so they can be multiplied for each likelihood evaluated in yij, j fixed
      for(int l=0; l<L; l++)
      { // we are creating 2 matrices where we have the subsample yj by row evaluated in theta_lk (different for each row). I need to sum over the different ms (so sum by column and the do the product)  
        Lik_forall_L.row(l)   = Likelihood_row(wrap(yj), tslk(arma::span(l), arma::span(k), arma::span::all));
        Lik_forall_L0.row(l)  = Likelihood_row(wrap(yj), tsl0.row(l));
      }
      
      WtimesLik  = Lik_forall_L % OSL;
      WtimesLik0 = Lik_forall_L0 % OSL0;
      
      arma::rowvec sommalk  = sum(WtimesLik,0),       
                   sommal0  = sum(WtimesLik0,0),    //sommo davvero per colonna?? si, sum(M,0) somma per colonna
                   sommal02 = sommal0*rho(k),
                   sommalk2 = sommalk*(1-rho(k));
      
      prob1[k] = sum( log( sommal02 + sommalk2 ));
    }

    prob2 = prob1 + log(psk);
    
    if(arma::is_finite(max(prob2))){
      prob3           = exp(prob2 - max(prob2));
      prob4           = prob3 / sum(prob3);
      actual_label[j] = RcppArmadillo::sample(possible_label, 1, 1, prob4)[0];  
    }else{
      Rcout << "All probabilities in subgroup " << j << " are zero :( \n";
      actual_label[j] = RcppArmadillo::sample(possible_label, 1, TRUE)[0];    
    }  }
  return  actual_label;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat  UPD_OBSER_LAB_XI(arma::colvec y_obser, 
                            arma::colvec y_group,
                            arma::cube tslk,   
                            arma::mat tsl0, 
                            arma::mat oslk,     
                            arma::colvec osl0, 
                            arma::colvec rho,
                            int J, int K, int L, int NN, 
                            arma::vec zj) {

  //data= 1st col: actual data, 2nd col actual group
    

  
  arma::colvec possible_label = arma::linspace<arma::vec>(1, 2*L, 2*L), 
               pp(2*L), pp2(2*L),  
    Chosen_Label(NN),Chosen_XI(NN); ;
  
    for(int j=0; j<J; j++){  //we do this for each actual group, in order to exploit vectorization where possible
    
    arma::uvec  inds = find( y_group == j+1 ) ;
    arma::vec     yj = y_obser.elem(inds);
    int n_yj         = size(yj)[0];
    arma::colvec obs_lab(n_yj), corresp_xi(n_yj);
    
    arma::mat 
    Lik_yij(2*L,n_yj), 
    OSL(L, n_yj), 
    OSL0(L, n_yj), 
    OSL_OSL0(2*L, n_yj),
    RHO(L, n_yj), 
    UmRHO(L, n_yj), 
    UmRHO_RHO(2*L, n_yj);
    
    RHO.fill(rho[zj[j]-1]); 
    UmRHO.fill(1-rho[zj[j]-1]);
    UmRHO_RHO = arma::join_cols(UmRHO,RHO); // 1-rho sta sopra, a moltiplicare la matrice di w_lk
    
    for(int ll=0; ll<(2*L); ll++){
        /*   attento a indici : prendi posizione giusta in zj con j, ma poi va riscalato indietro! */
      if(ll<L){
        Lik_yij.row(ll) = Likelihood_row(wrap(yj), tslk(arma::span(ll), arma::span( zj[j]-1 ), arma::span::all));
      }else{
        Lik_yij.row(ll) = Likelihood_row(wrap(yj), tsl0.row(ll-L));
      }
    }
    
    OSL  = arma::repelem(oslk.col(zj[j]-1), 1, n_yj);  
    OSL0 = arma::repelem(osl0, 1, n_yj);
    OSL_OSL0 = arma::join_cols(OSL,OSL0);
    arma::mat prob1 = log(Lik_yij) + log(OSL_OSL0) + log(UmRHO_RHO);
  
  // we obtain a 2L*nj matrix, where for each column we have the 2L probabilities corresponden to the various label Mij that we can assign to yij (here j fixed, to exploit vectorization) 
    
    for(int u=0; u<n_yj; u++){ //for each column we sample the observational label and the correspondent xi
      
      pp = prob1.col(u);
      
      if(arma::is_finite(max(pp))){
        pp2 = exp(pp - max(pp));
        obs_lab(u) = RcppArmadillo::sample(possible_label, 1, TRUE, pp2)[0];    
        if(obs_lab(u) >= (L+1)){ 
          corresp_xi(u)=1;
        }else{
          corresp_xi(u)=0;
        }      
      }else{
        Rcout << "All probabilities in subgroup " << j << "are zero :( \n";
        obs_lab(u) = RcppArmadillo::sample(possible_label, 1, TRUE)[0];    
        if(obs_lab(u) >= (L+1)){ 
          corresp_xi(u)=1;
        }else{
          corresp_xi(u)=0;
        }      
        }
    }
    Chosen_Label(inds) = obs_lab;
    Chosen_XI(inds)    = corresp_xi;
    }
  
  arma::mat Chosen_LabelXI = arma::join_rows(Chosen_Label,Chosen_XI);
  return Chosen_LabelXI;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::colvec my_rep_each(NumericVector x, NumericVector y){
  // given a vector of times, it repeats each element of the first vector
  int ni = x.size();
  NumericVector a;
  arma::colvec b, c;
  for(int i=0; i<ni; i++){
    a = rep(x[i],y[i]);
    c = as<arma::colvec>(a);
    b = arma::join_cols(b,c);
  }
  return(b);
}


// [[Rcpp::export]]
List Stratified_operations_lk(arma::vec x, 
                              arma::vec col1, int val1, 
                              arma::vec col2, int val2){
  
  arma::uvec inds1      = find( col1 == val1 ) ;
  arma::uvec inds2      = find( col2 == val2 ) ;
  arma::uvec ind12      = intersect(inds1, inds2);
  arma::vec  x12        = x.elem(ind12);
  
  int    n12  = x12.size(); 
  double xm12 = 0.0, 
         sm12 = 0.0;
  
  if(n12>0){   
    xm12 = mean(x12);
    sm12 = sum((x12-xm12)%(x12-xm12));
  }
  
  return List::create(
    _["ind12"] = ind12,
    _["x12"]   = x12,
    _["xm12"]  = xm12,
    _["sm12"]  = sm12,
    _["n12"]   = n12);
}


// [[Rcpp::export]]
List Stratified_operations_0(arma::vec x, 
                             arma::vec col1, int val1){
  
  arma::uvec inds1      = find( col1 == val1 ) ;
  arma::vec  x10        = x.elem(inds1);
  
  int n10 = x10.size(); 
  double xm10=0.0, sm10=0.0;
  if(n10 > 0){   
    xm10 = mean(x10);
    sm10 = sum((x10-xm10)%(x10-xm10));
  }
  
  return List::create(
    _["inds1"] = inds1,
    _["x10"]   = x10,
    _["xm10"]  = xm10,
    _["sm10"]  = sm10,
    _["n10"]   = n10);
}



// [[Rcpp::export]]
List Int_cpp(arma::colvec y_obser,
             arma::colvec y_group,
             arma::vec  cij, 
             NumericVector zj, 
             NumericVector nj, 
             int K, int L){
  
  arma::mat data = arma::join_rows(y_obser,y_group);
  List R;
  arma::vec lab_g = my_rep_each(zj,nj);
  arma::mat    M_m(L,K), M_s(L,K), A(L,K), nlk(L,K);
  data = arma::join_rows(arma::join_rows(data,cij),lab_g);
  
  for(int l =0; l<L; l++){
    for(int k =0; k<K; k++){
      R = Stratified_operations_lk(y_obser, lab_g, k+1, cij, l+1);    
      M_m(l,k) = R["xm12"];
      nlk(l,k) = R["n12"];
      M_s(l,k) = R["sm12"];
    } 
  }
  arma::colvec M_m0(L),  M_s0(L),  A0(L),  nl0(L);
  
  for(int ls = 0; ls<L; ls++){
      R = Stratified_operations_0(y_obser, cij, ls+1+L);    
      M_m0(ls) = R["xm10"];
      nl0(ls)  = R["n10"];
      M_s0(ls) = R["sm10"];
  }
  
  return(List::create(
      _["A"]   =  data,
      _["m"]   =  M_m,
      _["s"]   =  M_s,
      _["nlk"] =  nlk,
      _["m0"]   =  M_m0,
      _["s0"]   =  M_s0,
      _["nl0"] =  nl0));
}


// Stick Breaking simulator
// [[Rcpp::export]]
List SB(int N, double a, double b) {
 
  NumericVector V(N), pi(N);
  double tmV;
  
  V = rbeta(N,a,b);
  V[N-1]=1; pi[0]=V[0];
  NumericVector mV = 1-V;
  
  
  for(int h=1; h<N; h++){
    for(int l=0; l<(h); l++){ tmV *= mV[l]; }
    pi[h] = V[h]*tmV;
    tmV=1;
  }  
  
  return List::create(
    _["V"]  = V,
    _["mV"] = mV,
    _["pi"] = pi);
  
}


// SB simulator given the vector of Beta realizations 
List SB_given_u(NumericVector V) {
  
  double tmV;       int N = V.size();
  NumericVector pi(N);
  
  pi[0]=V[0];
  NumericVector mV = 1-V;
  
  for(int h=1; h<N; h++){
    tmV=1;
    for(int l=0; l<h; l++){ tmV *= mV[l]; }
    pi[h] = V[h]*tmV;
  }  
  
  return List::create(
    _["V"]  = V,
    _["mV"] = mV,
    _["pi"] = pi);
}



// A function that creates a cube useful for the posteriors of the SB weights
// [[Rcpp::export]]
arma::cube gimme_B(arma::mat  nlk, double beta) {
  
  int L = size(nlk)[0], K = size(nlk)[1];
  arma::mat A(L,L);
  
  for(int i=0; i<L; i++){
    for(int j=0; j<L; j++){                     //crea matrice di 1 tranne che sulla diagonale
      if(i==j){A(i,j)=0;}else{A(i,j)=1;}
    }
  }
  
  arma::mat At = trimatl(A),                    // estrae diagonale inferiore da matrice A --> At è matrice diagonale inferiore d 1
            NLK = nlk+1,
            seconda =  trans( trans(nlk) * At ) + beta  ;       
                                                // il trucco delle diagonale serve a sommare tutto ciò che appare da un certo indice in poi, utile per trovare gli argomenti della beta per le componenti dello stick breaking a posteriori 
  arma::cube B(L-1, K, 2);
  
  B.slice(0) = NLK.submat( 0, 0, L-2, K-1 );
  B.slice(1) = seconda.submat( 0, 0, L-2, K-1 ) ;
  
  return B;
}


// function that actually computes the SB weights a posteriori
// [[Rcpp::export]]
List updt_omegastar_lk(arma::mat nlk, double beta){
  int L = size(nlk)[0], K = size(nlk)[1];
  double a,b;
  
  arma::cube B = gimme_B(nlk, beta);
  arma::mat  v(L,K), os(L,K);
  List SSS;
  
  for(int g=0; g<L; g++){
    if(g==L-1){
      for(int j=0; j<K; j++){
        v(g,j) = 1;   
      }   
    }else{
      for(int j=0; j<K; j++){
         a = B(arma::span(g),arma::span(j),arma::span(0))[0];
         b = B(arma::span(g),arma::span(j),arma::span(1))[0];
        v(g,j) = R::rbeta(a,b);   
      }}}
  
  for(int t=0; t<K; t++){
    SSS = SB_given_u(wrap(v.col(t))); 
    os.col(t) = as<arma::colvec>(SSS["pi"]);
  }
  
  return List::create(
    _["os"]  = os,
    _["v"] = v);
}

// creates a table of a vector, keeping also the absences (zeros)
// [[Rcpp::export]]
arma::vec table_factor(arma::colvec X, arma::colvec levels){
  int n = levels.size(), nX = X.size();
  arma::vec tabbozzo(n),  count(nX);
  for(int i=0; i < n; i++){
    for(int j=0; j < nX; j++){
      count[j] =  (X[j]==levels[i]);
    }
    tabbozzo[i]=accu(count);
  } 
  return(tabbozzo);
}

// update the distributional stick breaking weights
// [[Rcpp::export]]
List updt_pistar_rcpp(arma::vec zj, int K, double alpha){
  
  arma::colvec lev = arma::linspace<arma::vec>(1, K, K),
               uk(K), os;
  arma::vec mk = table_factor(zj, lev);
  List SSS;
  
  for(int l=0; l<K; l++){
    if(l==(K-1)){
      uk[l] = 1;
    }else{
      arma::uvec  pos = arma::linspace<arma::uvec>(l+1, K-1, K-1-l);
      double a=1+mk[l],
             b= alpha + accu( mk.elem(pos) );
      uk[l] = R::rbeta(a,b); 
    }}
  
    SSS = SB_given_u(wrap(uk)); 
    os = as<arma::colvec>(SSS["pi"]);
    
    return List::create(
      _["pistar"]  = os,
      _["uk"] = uk);
    }


// update the distributional stick breaking weights
// [[Rcpp::export]]
List updt_omegastar_0(arma::vec nl0, int L, double alpha){
  
  arma::colvec uk(L), o0;
  List SSS;
  
  for(int l=0; l<L; l++){
    if(l==(L-1)){
      uk[l] = 1;
    }else{
      arma::uvec  pos = arma::linspace<arma::uvec>(l+1, L-1, L-1-l);
      double a=1+nl0[l],
             b= alpha + accu( nl0.elem(pos) );
      uk[l] = R::rbeta(a,b); 
    }}
  
    SSS = SB_given_u(wrap(uk)); 
    o0 = as<arma::colvec>(SSS["pi"]);
    
    return List::create(
      _["pistar"]  = o0,
      _["uk"] = uk);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec upd_rho_j(arma::vec xi, 
                    arma::colvec lab_g, 
                    double A, double B, int K){
  
  arma::colvec newRho(K);
  
  for(int k = 0; k < K; k++){
    arma::uvec ind = find( lab_g==(k+1) ); // qui vi era buggone!
    int n_k        = ind.n_rows;
    int som        = sum(xi.elem(ind));
    newRho[k]      = R::rbeta(som + A, n_k + B - som);
  }
  return(newRho);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double rt_cpp(double nu, double lambda){
  double TAA;
  TAA = R::rnorm(lambda,1.0) / exp(.5*log( R::rchisq(nu)/nu ));
  return(TAA);
}


// [[Rcpp::export]]
arma::cube UPD_TH_CD_cpp_T_lk( List AIFD,  
                               double k0, double m0, 
                               double b0, double a0, 
                               int L, int K){
  
  arma::cube tslk(L,K,2); tslk.fill(0);
  
  arma::mat AIFD_m = as<arma::mat>(AIFD["m"]);
  arma::mat AIFD_s = as<arma::mat>(AIFD["s"]);
  arma::mat nlk    = as<arma::mat>(AIFD["nlk"]);
  
  arma::mat mu_star_star = (k0*m0+AIFD_m%nlk)%exp(-1*log(k0+nlk));
  arma::mat k_star_star  = k0 + nlk;
  arma::mat a_star_star  = a0 + nlk/2;
  arma::mat b_star_star  = b0 + .5*AIFD_s + ((k0 * nlk) % pow( (AIFD_m-m0),2)) % pow((2*(k0+nlk)),-1);

    for( int l=0; l<L; l++){
    for( int k=0; k<K; k++){ 
      
      double isigma = 1/ rgamma(1,a_star_star(l,k),1/(b_star_star(l,k)))[0];
      tslk(arma::span(l), arma::span(k), arma::span(1)) = isigma;
      tslk(arma::span(l), arma::span(k), arma::span(0)) = 
       rt_cpp(2*a_star_star(l,k), 0.0) *
       pow( ( b_star_star(l,k)/(k_star_star(l,k)*a_star_star(l,k)) ),.5)+
       mu_star_star(l,k);
    }
  }  
  return(tslk); 
}

// [[Rcpp::export]]
arma::mat UPD_TH_CD_cpp_T_0( List AIFD,  
                             double k0, double m0, 
                             double a0, double b0, 
                             int L, int K){
  
  arma::mat tsl0(L,2); tsl0.fill(0);
  
  arma::mat AIFD_m0 = as<arma::mat>(AIFD["m0"]);
  arma::mat AIFD_s0 = as<arma::mat>(AIFD["s0"]);
  arma::mat nl0     = as<arma::mat>(AIFD["nl0"]);
  
  arma::mat mu_star_star = (k0*m0+AIFD_m0%nl0)%exp(-1*log(k0+nl0));
  arma::mat k_star_star  = k0 + nl0;
  arma::mat a_star_star  = a0 + nl0/2;
  arma::mat b_star_star  = b0 + .5 * AIFD_s0 + (k0 * nl0 % pow( (AIFD_m0 - m0),2) ) % pow( ( 2 * (k0+nl0) ), -1);

    for( int l=0; l<L; l++){

      double isigma = 1/ rgamma(1,a_star_star(l),1/(b_star_star(l)))[0];
      tsl0(arma::span(l), arma::span(1)) = isigma;
      tsl0(arma::span(l), arma::span(0)) = 
       rt_cpp(2*a_star_star(l), 0.0) *
       pow( ( b_star_star(l)/(k_star_star(l) * a_star_star(l)) ),.5)+
       mu_star_star(l);
    }
  return(tsl0); 
}


// [[Rcpp::export]]
arma::cube heat_mat(arma::mat NANA) {
  int  J = NANA.n_rows , 
    NSIM = NANA.n_cols;
  arma::cube CU(J, J, NSIM);
  
  for(int i=0; i<NSIM ;i++){
    for(int l=0; l<J ;l++){
      for(int k=0; k<J ;k++){
        if( NANA(l,i)==NANA(k,i) ){
          CU(l, k, i) = 1;
        }else{
          CU(l, k, i) = 0; } 
      }}}
  return CU;
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
List Ciclone(int NSIM, int burnIN,
             arma::colvec y_obser, arma::colvec y_group, 
             arma::mat data,  
             int  K, int L, int J, int N,
             arma::cube tslk, arma::mat tsl0, 
             arma::colvec psk, arma::vec rho,
             arma::cube oslk, arma::mat osl0,
             double alpha, double beta){

  arma::vec zj(J), xi(N), cij(N), uk(K);
  arma::mat MAMA(N,2);
  List res_pistar;

  
  
  for(int sim=0; sim<NSIM; sim++){ 
  
  
    zj = UPD_DISTR_LAB(y_obser = y_obser,y_group = y_group, tslk = tslk, tsl0 = tsl0, psk = psk, rho = rho, oslk = oslk, osl0 = osl0, J =  J, K = K, L = L);
    ///////////////////////////////////////////////////////////////////////////////////////////////
    MAMA = UPD_OBSER_LAB_XI(y_obser = y_obser,y_group = y_group, tslk = tslk, tsl0 = tsl0, oslk = oslk, osl0 = osl0, rho = rho, J =J, L=L, K=K, zj=zj);
    xi   = MAMA.col(1);
    cij  = MAMA.col(0);
    //////////////////////////////////////////////////////////////////////////////////////////////////
    res_pistar  = updt_pistar_rcpp(zj=zj,alpha = alpha,K = K);
      psk       = res_pistar["pistar"];
      uk        = res_pistar["uk"];
  
  }
      
####################################################################################################################################################################
# INTERVAL: updating useful quantities
####################################################################################################################################################################
      
      AIFD <- Int_cpp(data = data.matrix(data), cij = cij, zj=zj, nj=nj, K=K, L=L)
        
####################################################################################################################################################################
# STEP 4, uptading w^*_lk
####################################################################################################################################################################
        
        res_ostar <- updt_omegastar_lk(nlk=AIFD$nlk,beta = beta)
          oslk      <- res_ostar$os
          v         <- res_ostar$v
          
          res_ostar0 <- updt_omegastar_0(nl0=AIFD$nl0,L = L, alpha = beta)
          osl0       <- res_ostar0$pistar
          v0         <- res_ostar0$uk
          
####################################################################################################################################################################
# Step 5  --  Update Theta
####################################################################################################################################################################
          
          tslk <- UPD_TH_CD_cpp_T_lk(AIFD = AIFD,k0 = k0, m0 = m0, b0 = b0, a0 = a0, L = L,K = K)
            tsl0 <- UPD_TH_CD_cpp_T_0( AIFD = AIFD,k0 = k0, m0 = m0, b0 = b0, a0 = a0, L = L,K = K)
            
####################################################################################################################################################################
# Step 6 - Simply updating alpha and beta
####################################################################################################################################################################
            
            if(fixedAB & sim == 2){
              alpha <- prior$alpha_DP ; 
              beta  <- prior$beta_DP
            }else{
              alpha <- 1#rgamma(1, a_alpha+(K-1),   b_alpha -  sum(log(1-uk[-K])))
              beta  <- 1#rgamma(1, a_beta +K*(L-1), b_beta  -  sum(log(1-v[-L,])))
            }
            
####################################################################################################################################################################
# Step 7 - Update Rhos
####################################################################################################################################################################
            
            if(RDG){
              rho <- rep(0,K)
            }else{
              rho <- upd_rho_j(xi = xi,lab_g = AIFD$A[,4], A = rho_A, B = rho_B, K = K)
            }
            
################################################################################################################################
# Storing the results
################################################################################################################################
            
            if (sim > burn_in & ((sim - burn_in) %% thinning == 0)) {
              rr                      <- floor((sim - burn_in)/thinning)
              pi_star_k[,rr]          <- psk
              omega_star_lk[,,rr]     <- oslk
              theta_star_lk[,,,rr]    <- tslk
              Csi_ij[,rr]             <- cij
              Z_j[,rr]                <- zj
              Xi[rr,]                 <- xi
              Omega_l0[,rr]           <- osl0
              theta_star_zero[,,rr]   <- tsl0
              Rho[rr,]                <- rho
#  if(sim %% m.step==0) plot(data$y_obser, col=rep(zj,nj)) #plot(BETAREG,type="l")
            }
# Printing the results
            if(sim==NSIM){
              RL <- list(psk  =  psk,
                         oslk =  oslk,
                         tslk =  tslk,
                         cij  =  cij,
                         zj   =  zj,
                         xi=xi, osl0=osl0, rho=rho )
            }
            
            if (verbose) {
              if (sim%%(verbose.step*thinning) == 0) {
                cat(paste("Sampling iteration: ", sim, " out of ",NSIM*thinning + burn_in, "---", length(unique(zj)), "\n",
                          sep = "" ))}}
}




} // closing Ciclone
 
 
 */