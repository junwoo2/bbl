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

// [[Rcpp::export]]
NumericMatrix predict_class(IntegerMatrix xid, IntegerVector Ly, List h, List J,
                LogicalVector numericmodel, NumericVector lz, NumericVector py){

  int n = xid.nrow();
  int m = xid.ncol();

  int ly = Ly[0];
  NumericMatrix ay(n,ly);
  
  for(int k=0; k<n; k++){
    std::vector<double> E(ly);
    for(int iy=0; iy<ly; iy++){
      double e=0;
      NumericMatrix hy = h[iy];
      List Jy = J[iy];
      for(int i=0; i<m; i++){
        if(xid(k,i)==0) continue;
        if(numericmodel[0]) 
          e += hy(i,0)*xid(k,i);
        else 
          e += hy(i,xid(k,i)-1);
        List Ji = Jy[i];
        for(int j=0; j<m; j++){
          if(j==i || xid(k,j)==0) continue;
          NumericMatrix Jj = Ji[j];
          if(numericmodel[0]) e += Jj(0,0)*xid(k,i)*xid(k,j)/2.0;
          else e += Jj(xid(k,i)-1,xid(k,j)-1)/2.0;
        }
      }
      E[iy] = e - lz[iy] + log(py[iy]);
    }
    for(int iy=0; iy<ly; iy++){
      double sum=0.0;
      for(int y2=0; y2<ly; y2++)
        if(y2!=iy) sum += exp(E[y2] - E[iy]);
      ay(k,iy) = -log(sum);
    }
  }

  return ay;
}
