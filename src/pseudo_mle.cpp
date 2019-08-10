#include <Rcpp.h>
#include <gsl/gsl_vector.h>
#include "bfgs.h"
using namespace Rcpp;

// [[Rcpp::export]]
List pseudo_mle(NumericMatrix xi, IntegerVector L, LogicalVector Numeric,
                NumericVector Lambda, IntegerVector Nprint, IntegerVector Itmax, 
                NumericVector Tol, IntegerVector Verbose){
  
  int n = xi.nrow();
  int m = xi.ncol();
  std::vector<std::vector<short> > sv(n);
  for(int i=0; i<n; i++) for(int j=0; j<m; j++)
    sv[i].push_back(xi(i,j));
  
  std::vector<std::vector<double> > h(m);
  std::vector<std::vector<std::vector<double> > > J(m);
  
  int Lv = L[0]-1;
  double lambda = Lambda[0];
  int nprint = Nprint[0];
  unsigned int Imax = Itmax[0];
  double tol = Tol[0];
  int verbose = Verbose[0];
  bool numeric = Numeric[0];
  
  double lkl = 0;
  double lz = 0;
  for(int i0=0; i0<m; i0++){
    double z=0;
    lkl += lpr_psl(i0, sv, Lv, lambda, h[i0], J[i0], nprint, Imax, tol,
                   verbose, z, numeric);
    lz += z;
  }
  lkl /= n;
  
  List x = List::create(Named("h") = h, Named("J") = J, 
                        Named("lkl") = lkl, Named("lz") = lz);
  return x;
}
