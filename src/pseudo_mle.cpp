#include <Rcpp.h>
#include <gsl/gsl_vector.h>
#include "bfgs.h"
using namespace Rcpp;

// [[Rcpp::export]]
List pseudo_mle(NumericMatrix xi, LogicalVector Numeric,
                NumericVector Lambda, IntegerVector Nprint, 
                IntegerVector Itmax, NumericVector Tol,
                LogicalVector Naive, IntegerVector Verbose){
  
  int n = xi.nrow();
  int m = xi.ncol();
  std::vector<std::vector<short> > sv(n);
  for(int i=0; i<n; i++) for(int j=0; j<m; j++)
    sv[i].push_back(xi(i,j));
  std::vector<short> L(m);
  std::vector<bool> bad(m);
  for(int i=0; i<m; i++){
    short xmin=0;
    short xmax=0;
    for(int k=0; k<n; k++){ 
      if(xmax<xi(k,i)) xmax=xi(k,i);
      if(xmin>xi(k,i)) xmin=xi(k,i);
    }
    L[i]=xmax;
    if(xmax==xmin) bad[i]=true;
    else bad[i]=false;
  }
  
  std::vector<std::vector<double> > h(m);
  std::vector<std::vector<std::vector<double> > > J(m);
  
  double lambda = Lambda[0];
  int nprint = Nprint[0];
  unsigned int Imax = Itmax[0];
  double tol = Tol[0];
  int verbose = Verbose[0];
  bool numeric = Numeric[0];
  bool naive = Naive[0];
  
  double lkl = 0;
  double lz = 0;
  for(int i0=0; i0<m; i0++){
    double z=0;
    bool failed=false;
    if(!bad[i0]){
      lkl += lpr_psl(i0, sv, L, lambda, h[i0], J[i0], nprint, Imax, tol,
                   verbose, z, numeric, naive, failed);
      if(failed)
        std::cout << " Warning: failed to converge in pseudo\n";
    }
    else{
      h[i0].push_back(0);
      J[i0].resize(m);
      for(int j=0;j<m;j++) J[i0][j].push_back(0);
    }
    lz += z;
    Rcpp::checkUserInterrupt();
  }
  lkl /= n;
  
  List x = List::create(Named("h") = h, Named("J") = J, 
                        Named("lkl") = lkl, Named("lz") = lz);
  return x;
}

// [[Rcpp::export]]
NumericVector predict_class(IntegerVector xid, IntegerVector Ly, List h, List J,
                LogicalVector numericmodel, NumericVector lz, NumericVector py,
                LogicalVector Naive){

  int m = xid.length();
  int ly = Ly[0];
  
  NumericVector E(ly);
  for(int iy=0; iy<ly; iy++){
    double e=0;
    List hy = h[iy];
    List Jy = J[iy];
    for(int i=0; i<m; i++){
      if(xid(i)==0) continue;
      NumericVector hi = hy[i];
      if(numericmodel[0]) 
        e += hi(0)*xid(i);
      else if(hi.length() < xid(i)) continue;
      else
        e += hi(xid(i)-1);
      List Ji = Jy[i];
      if(Naive[0]) continue;
      for(int j=0; j<m; j++){
        if(j==i || xid(j)==0) continue;
        NumericMatrix Jj = Ji[j];
        if(numericmodel[0]) e += Jj(0,0)*xid(i)*xid(j)/2.0;
        else if(Jj.nrow()<xid(i) ||
                Jj.ncol()<xid(j)) continue;
        else e += Jj(xid(i)-1,xid(j)-1)/2.0;
      }
    }
    E[iy] = e - lz[iy] + log(py[iy]);
    Rcpp::checkUserInterrupt();
  }
  
  return E;
}
