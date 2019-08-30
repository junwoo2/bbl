#include <cmath>
#include <Rcpp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include "bfgs.h"

using namespace std;

void invC(const vector<vector<short> > &ai, const vector<short> &L, 
          double &lnz, vector<vector<double> > &h,
          vector<vector<vector<double> > > &J, double eps){

  int nsnp=L.size();
  int ndim=0;
  vector<vector<double> > f1(nsnp);
  vector<vector<vector<double> > > f2(nsnp);
  for(int i=0; i<nsnp; i++){
    f12(i, ai, f1[i], f2[i], L, false, true);
    ndim += L[i];
  }

  gsl_matrix *A;
  gsl_matrix *Ai;
  gsl_permutation *perm;

  A=gsl_matrix_alloc(ndim,ndim);   // serial version using GSL
  Ai=gsl_matrix_alloc(ndim,ndim);
  perm=gsl_permutation_alloc(ndim);

  double tr=0;
  for(int i=0;i<nsnp;i++) for(int l=0;l<L[i];l++)
    tr+=f1[i][l]*(1-f1[i][l]);
  tr/=ndim;

  int idx=0;
  for(int i=0;i<nsnp;i++) for(int l0=0;l0<L[i];l0++){
    int jdx = 0;
    for(int j=0;j<nsnp;j++) for(int l1=0;l1<L[j];l1++){
      double x=eps*(f2[i][j][L[j]*l0+l1]-f1[i][l0]*f1[j][l1]);
      if(i==j && l0==l1)
        x += (1-eps)*tr;
      gsl_matrix_set(A, idx, jdx++, x);
    }
    idx++;
  }

  int s;
  gsl_linalg_LU_decomp(A,perm,&s);
  gsl_linalg_LU_invert(A,perm,Ai);

  h.resize(nsnp);
  J.resize(nsnp);
  lnz=0;
  idx=0;
  for(int i=0;i<nsnp;i++){
    h[i].resize(L[i]);
    J[i].resize(nsnp);
    for(int j=0;j<nsnp;j++) J[i][j].resize(L[i]*L[j]);
    double s0=0;
    for(int l0=0; l0<L[i]; l0++) s0 += f1[i][l0];
    lnz+=-log(1-s0);
    for(int l0=0;l0<L[i];l0++){
      double f=log(f1[i][l0]/(1.0-s0));
      int jdx=0;
      for(int j=0;j<nsnp;j++) for(int l1=0;l1<L[j];l1++){
        if(i!=j){
          double x=gsl_matrix_get(Ai, idx, jdx);
          J[i][j][L[j]*l0+l1]=-x;
          f+=x*f1[j][l1];
          lnz+=0.5*x*f1[i][l0]*f1[j][l1];
        }
        jdx++;
      }
      h[i][l0]=f;
      idx++;
    }
  }

  gsl_matrix_free(A);
  gsl_matrix_free(Ai);
  gsl_permutation_free(perm);

}
