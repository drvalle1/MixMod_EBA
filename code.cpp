#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
// code to sample from categorical distribution
int whichLessDVPresenceFast(double value, 
                            NumericVector prob) {
  int res=-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob[i];
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

// [[Rcpp::export]]
// This function samples z's
IntegerVector samplez(NumericMatrix ltheta,
                      IntegerMatrix nmat,
                      IntegerMatrix Nminusn,
                      NumericMatrix lpi1,
                      NumericMatrix l1minuspi1,
                      int ngroup,
                      int nloc,
                      int nheight,
                      NumericVector randu,
                      IntegerVector TransID) {

  IntegerVector zvec(nloc);
  NumericVector prob(ngroup);
  int znew;
  double soma=0;
   
  for(int i=0; i<nloc;i++){
    for (int k=0; k<ngroup; k++){
      soma=0;
      for (int j=0; j<nheight; j++){
        soma=soma+nmat(i,j)*lpi1(k,j)+
                  Nminusn(i,j)*l1minuspi1(k,j);
      }
      prob[k]=soma+ltheta(TransID[i],k);
    }
    prob=prob-max(prob);
    prob=exp(prob);
    prob=prob/sum(prob);

    //multinomial draw
    znew=whichLessDVPresenceFast(randu[i],prob);
    zvec[i]=znew+1;
  }
  return zvec;
}

// [[Rcpp::export]]
// This function calculates ntk
IntegerMatrix calc_ntk(IntegerVector z,
                      int ngroup,
                      int nloc,
                      int ntransect,
                      IntegerVector TransID) {
  
  IntegerMatrix ntk(ntransect,ngroup);

  for(int i=0; i<nloc;i++){
    ntk(TransID[i],z[i])=ntk(TransID[i],z[i])+1;
  }
  return ntk;
}