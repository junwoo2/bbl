#include <Rcpp.h>
#include <gsl/gsl_vector.h>
#include "bfgs.h"
using namespace Rcpp;

// [[Rcpp::export]]
List pseudo_mle(NumericMatrix xi, IntegerVector L, NumericVector Lambda){
  
  int n = xi.nrow();
  int m = xi.ncol();
  std::vector<std::vector<short> > sv(n);
  for(int i=0; i<n; i++) for(int j=0; j<m; j++)
    sv[i].push_back(xi(i,j));
  
  std::vector<std::vector<double> > h(m);
  std::vector<std::vector<std::vector<double> > > J(m);
  
  int Lv = L[0]-1;
  double lambda = Lambda[0];
  
  double lkl = 0;
  for(int i0=0; i0<m; i0++)
    lkl += lpr_psl(i0, sv, Lv, lambda, h[i0], J[i0]);
  lkl /= n;
  
  List x = List::create(Named("h") = h, Named("J") = J, 
                        Named("lkl") = lkl);
  return x;
}
