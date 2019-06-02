#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;



// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat distMat(int nn, int dd, arma::mat ZZ)
{
  arma::mat dMat(nn,nn,arma::fill::zeros);
  int ii,jj,kk;
  double tmp;
  for(ii = 0 ; ii <= (nn-1) ; ii++){
    for(jj = 0 ; jj <= ii ; jj++){
      tmp = 0.0;
      for(kk = 0 ; kk < dd ; kk++){
        //    double VV = ZZ(ii+kk*nn) - ZZ(jj+kk*nn);
        double VV = ZZ(ii,kk) - ZZ(jj,kk);
        tmp = tmp + VV*VV;
      }
      dMat(ii,jj) = sqrt(tmp);
      dMat(jj,ii) = sqrt(tmp);
    }
  }
  return(dMat);
}

double logitInverse(double x){
  return 1.0/(1.0 + exp(-1.0 * x));
}



// [[Rcpp::export]]
double FullLogLikCOVSR(arma::mat YY, arma::mat ZZ,arma::mat XX,
                     arma::vec Beta, arma::vec SS, arma::vec RR,
                     double intercept,int nn,int dd, int pp){
  arma::mat dMat = distMat(nn,dd,ZZ);
  double total = 0.0;
  double tmp1,tmp2;
  for(int ii = 1; ii < nn; ii++){ //we exclude diagonal elements
    for(int jj = 0; jj < ii; jj++){
      double tmpVal1 = 0.0;
      double tmpVal2 = 0.0; 	
      int nacountij = 0;
      int nacountji = 0;	
      int Ri = RR(ii);
      int Rj = RR(jj);
      int Si = SS(ii);
      int Sj = SS(jj);
      for(int kk = 0; kk < pp; kk++){
        double Xij = XX(ii+kk*nn,jj);
        double Xji = XX(jj+kk*nn,ii);
        if(R_IsNA(Xij)){
          nacountij = nacountij + 1;
        }else{
          tmpVal1 = tmpVal1 + Beta(kk)*Xij;
        } 
        if(R_IsNA(Xji)){nacountji = nacountji +1;
        }else{		
          tmpVal2 = tmpVal2 + Beta(kk)*Xji;
        }
      }
      if(nacountij == 0){	
        double vv1 = (intercept + tmpVal1 - dMat(ii,jj) + Si +Rj); 
        tmp1 = logitInverse(vv1); 
      }
      if(nacountji == 0){	
        double vv2 = (intercept + tmpVal2 - dMat(jj,ii)+ Sj+Ri);
        tmp2 = logitInverse(vv2);
      }
      if(R_IsNA(YY(ii,jj))|| nacountij > 0){
        total = total + 0.0;
      }else{
        if(YY(ii,jj) == 1){
          total = total + log(tmp1);
        }else if(YY(ii,jj) == 0){
          total = total + log(1.0 - tmp1);
        }       }			
      if(R_IsNA(YY(jj,ii))|| nacountji > 0){
        total = total + 0.0;
      }else{
        if(YY(jj,ii) == 1){
          total = total + log(tmp2);
        }else if(YY(jj,ii) == 0){
          total = total + log(1.0 - tmp2);
        }   
      }
    }		
  }
  return total;
}



// [[Rcpp::export]]
double likelihoodiCOVSR(int ii,int dd,int nn, int pp, arma::mat Yt,
                      arma::mat Xt,arma::mat Zt,
                      arma::vec SS, arma::vec RR,
                      double intercept,arma::vec Beta)
{
  ii = ii - 1.0;
  double lliki = 0.0;
  double pij,pji;
  arma::mat dMat =  distMat(nn,dd,Zt);
  int Si = SS(ii);
  int Ri = RR(ii);
  //  dMat.print("dMat:");
  //  #compute loglikelihood for the entries of Yt except the diagonal
  for(int jj =0; jj < nn; jj++){
    if(jj == ii){
      lliki = lliki + 0.0;
    }
    if(jj != ii){ 
      int Rj = RR(jj);
      int Sj = SS(jj);
      double tmpVal1 = 0.0;
      double tmpVal2 = 0.0;
      int nacountij = 0;
      int nacountji = 0;	
      for(int kk = 0; kk < pp; kk++){
        double Xij = Xt(ii+kk*nn,jj);
        double Xji = Xt(jj+kk*nn,ii);
        if(R_IsNA(Xij)){nacountij = nacountij + 1;
        }else{
          tmpVal1 = tmpVal1 + Beta(kk)*Xt(ii+kk*nn,jj);
        } 
        if(R_IsNA(Xji)){nacountji = nacountji +1;
        }else{
          tmpVal2 = tmpVal2 + Beta(kk)*Xt(jj+kk*nn,ii);
        }		}
      if(nacountij == 0){	
        double vv1 = (intercept + tmpVal1 - dMat(ii,jj)+Si+Rj); 
        pij = logitInverse(vv1); 
      }
      if(nacountji == 0){	
        double vv2 = (intercept + tmpVal2 - dMat(jj,ii)+Sj+Ri);
        pji = logitInverse(vv2);
      }
      if(R_IsNA(Yt(ii,jj))|| nacountij > 0){
        lliki = lliki + 0.0;
      }else{	
        if(Yt(ii,jj) == 1){
          lliki = lliki + (Yt(ii,jj))*log(pij);
        }
        if(Yt(ii,jj) == 0){
          lliki = lliki + (1-Yt(ii,jj))*log(1-pij);
        }
      }
      if(R_IsNA(Yt(jj,ii))|| nacountji > 0){
        lliki = lliki + 0.0;
      }else{	
        if(Yt(jj,ii) == 1){
          lliki = lliki + (Yt(jj,ii)*log(pji));
        }
        if(Yt(jj,ii) == 0){
          lliki = lliki + (1-Yt(jj,ii))*log(1-pji);
        }   
      }
    } }
  return lliki;
}    


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

