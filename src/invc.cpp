#include <cmath>
#include <Rcpp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include "bfgs.h"

const int Lbig=10;

using namespace std;

double fnumeric(int L, double f){
  
  double x = log((L-f)/(L-f-1.0));
//  double x=0.1;
  int maxit=100;
  double tol=1.0e-3;
  
  if(fabs(f-0.5)<tol) return(0.0);
  
  if(L>=Lbig && f>0) return(x);
  int i=0;
  for(i=0; i<maxit; i++){
    double elx = exp(-L*x);
    double ex = exp(-x);
    double fum = L/(1.0-elx) - 1.0/(1.0-ex) - f;
    double dfum = -L*L*elx/(1.0-elx)/(1.0-elx) + ex/(1.0-ex)/(1.0-ex);
//    double fum = f*elx*ex-(1+f)*elx+(L-f)*ex+f-L+1;
//    double dfum = -f*(L+1)*elx*ex+L*(1+f)*elx-(L-f)*ex;
    double xp = x - fum/dfum;
    double df = fabs((xp/x)*(xp/x)-1.0);
    if(df < tol) break;
    x = xp;
  }
  if(i >= maxit){
    Rcpp::Rcerr << "Maximum iteration limit reached\n";
    Rcpp::Rcerr << x <<" " << f << endl;
  }

  return(x);
}
void invC(const vector<vector<short> > &ai, const vector<int> &frq,
          const vector<bool> &numeric, const vector<short> &L, 
          double &E, double &lnz, vector<vector<double> > &h,
          vector<vector<vector<double> > > &J, double eps){

  int nsnp=L.size();
  int ndim=0;
  vector<vector<double> > f1(nsnp);
  vector<vector<vector<double> > > f2(nsnp);
  for(int i=0; i<nsnp; i++){
    f12(i, ai, frq, numeric, f1[i], f2[i], L, false, true);
    if(numeric[i]) ndim++;
    else ndim += L[i];
  }
//  cout << f2[0][0][0] << endl;

  gsl_matrix *A;
  gsl_matrix *Ai;
  gsl_permutation *perm;

  if(eps>0){
    A=gsl_matrix_alloc(ndim,ndim);   // serial version using GSL
    Ai=gsl_matrix_alloc(ndim,ndim);
    perm=gsl_permutation_alloc(ndim);
  }

  double tr=0;
  for(int i=0;i<nsnp;i++){ 
    int Li=L[i];
    if(numeric[i]) Li=1;
    for(int l=0;l<Li;l++){
      if(numeric[i])
        tr += f2[i][i][0]-f1[i][0]*f1[i][0];
      else
        tr+=f1[i][l]*(1-f1[i][l]);
    }
  }
  tr/=ndim;
//  cout << tr << endl;

  int idx=0;
  for(int i=0;i<nsnp;i++){
    int Li= (numeric[i]? 1 : L[i]);
    for(int l0=0;l0<Li;l0++){
      int jdx = 0;
      for(int j=0;j<nsnp;j++){ 
        int Lj = (numeric[j]? 1 : L[j]);
        for(int l1=0;l1<Lj;l1++){
          double x=eps*(f2[i][j][Lj*l0+l1]-f1[i][l0]*f1[j][l1]);
          if(eps==0.0) continue;
          if(i==j && l0==l1)
            x += (1-eps)*tr;
          gsl_matrix_set(A, idx, jdx++, x);
//          cout << x << endl;
        }
      }
      idx++;
    }
  }

  int s;
  if(eps>0){
    gsl_linalg_LU_decomp(A,perm,&s);
    gsl_linalg_LU_invert(A,perm,Ai);
  }

  h.resize(nsnp);
  J.resize(nsnp);
  lnz=0;
  idx=0;
  for(int i=0;i<nsnp;i++){
    int Li= (numeric[i]? 1 : L[i]);
    h[i].resize(Li);
    J[i].resize(nsnp);
    for(int j=0;j<nsnp;j++){
      int Lj = (numeric[j]? 1 : L[j]);
      J[i][j].resize(Li*Lj);
    }
    double f=0.0;
    double s0=0.0;
    if(numeric[i]){
//      cout << f1[i][0] << endl;
      f = fnumeric(L[i]+1, f1[i][0]);
      if(L[i]+1 < Lbig){
        if(fabs(f)<Tol) lnz += L[i]+1;
        else lnz += log((1.0-exp((L[i]+1)*f))/(1.0-exp(f)));
      }
      else if(f>0)
        lnz += (L[i]+1)*f - log(exp(f)-1.0);
      else
        lnz += log(fabs((L[i]+1)*f)/(1.0-exp(f)));
    }
    else{
      for(int l0=0; l0<Li; l0++) s0 += f1[i][l0];
      lnz += -log(1-s0);
    }
    for(int l0=0;l0<Li;l0++){
      if(!numeric[i]) 
        f=log(f1[i][l0]/(1.0-s0));
//      cout << f1[i][0] << endl;
      if(eps>0.0){
        int jdx=0;
        for(int j=0;j<nsnp;j++){
          int Lj = (numeric[j]? 1 : L[j]);
          for(int l1=0;l1<Lj;l1++){
            if(i!=j){
              double x=gsl_matrix_get(Ai, idx, jdx);
              J[i][j][Lj*l0+l1]=-x;
              f+=x*f1[j][l1];
              lnz += 0.5*x*f1[i][l0]*f1[j][l1];
            }
          }
          jdx++;
        }
      }
      h[i][l0]=f;
      idx++;
    }
  }
  E=0;
  int nind=ai.size();
  for(int k=0; k<nind; k++){
    for(int i=0; i<nsnp; i++){
      int a0=ai[k][i];
      if(a0==0) continue;
      if(numeric[i])
        E += h[i][0];
      else
        E += h[i][a0-1];
      for(int j=i+1;j<nsnp;j++){
        int a1=ai[k][j];
        if(a1==0) continue;
        if(numeric[i] & numeric[j])
          E += J[i][j][0];
        else if(numeric[i])
          E += J[i][j][L[j]*(a0-1)];
        else if(numeric[j])
          E += J[i][j][a1-1];
        else
          E += J[i][j][L[j]*(a0-1)+a1-1];
      }
    }
  }
  E = E/nind - lnz;

  if(eps>0){
    gsl_matrix_free(A);
    gsl_matrix_free(Ai);
    gsl_permutation_free(perm);
  }
}
