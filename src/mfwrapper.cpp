#include <Rcpp.h>
#include <gsl/gsl_vector.h>
#include "bfgs.h"
using namespace Rcpp;

// [[Rcpp::export]]
List mfwrapper(NumericMatrix xi, IntegerVector Lv, NumericVector Eps){
  
  int n = xi.nrow();
  int m = xi.ncol();
  std::vector<short> L;
  std::vector<std::vector<short> > sv(n);
  for(int i=0; i<m; i++){
    L.push_back(Lv[i]);
    for(int k=0; k<n; k++)
      sv[k].push_back(xi(k,i));
  }
  
  int mL=L.size();
  std::vector<std::vector<double> > hp(mL);
  std::vector<std::vector<std::vector<double> > > Jp(mL);
  
  double eps = Eps[0];
  
  double lnz = 0;
  invC(sv, L, lnz, hp, Jp, eps);
  
  std::vector<std::vector<double> > h(m);
  std::vector<std::vector<std::vector<double> > > J(m);
  
  int ic=0;
  for(int i=0; i<m; i++){
    J[i].resize(m);
    h[i]=hp[ic];
    int jc=0;
    for(int j=0; j<m; j++)
      J[i][j]=Jp[ic][jc++];
    ic++;
  }
    
  List x = List::create(Named("h") = h, Named("J") = J, Named("lz") = lnz);
  return x;
}
