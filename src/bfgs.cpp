#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include "bfgs.h"

using namespace std;

void pan3(vector<double> &peff, int nsnp, int i0, int L, 
          const vector<short> &ci, vector<double> h1, 
          const vector<vector<double> > &J1){

  peff.resize(L);
  for(int a=0;a<L;a++){
    double e=h1[a];
    for(int j=0;j< nsnp;j++){
      if(j==i0) continue;
      int b=ci[j];
      if(b==0) continue;
      e+=J1[j][L*a+b-1];
    }
    peff[a]=e;
  }
  double max=0;
  for(int a=0;a<L;a++)
    if(peff[a]>max) max=peff[a];
  double z = exp(-max);
  for(int a=0;a<L;a++){
    peff[a] = exp(peff[a]-max);
    z+=peff[a];
  }
  for(int a=0;a<L;a++)
    peff[a]/=z;

}

double pan2(int nsnp, int i0, int L, const vector<short> &ci, 
  const vector<double> &h1, const vector<vector<double> > &J1){

  vector<double> peff(L);
  pan3(peff, nsnp, i0, L, ci, h1, J1);

  int a0=ci[i0];
  if(a0>0)
    return peff[a0-1];
  double p=1;
  for(int l=0; l<L; l++)
    p -= peff[l];
   
  return p;
}

double lnl_psl(const gsl_vector *v,void *params){  // evaluates log likelihood

  double ln;
  Param *par=(Param *)params;
  int i0=par->i0;
  int L=par->L;
  double lambda=par->lambda;
  int nsnp=(par->ai)[0].size();

  vector<double> h1(L);
  vector<vector<double> > J1(nsnp);
  for(int i=0; i<nsnp; i++) J1[i].resize(L*L);

  int m=0;
  for(int l0=0;l0<L;l0++){
    h1[l0]=gsl_vector_get(v,m++);
    for(int i=0; i<nsnp; i++){ 
      if(i==i0) continue;
      for(int l1=0;l1<L;l1++)
        J1[i][L*l0+l1]=gsl_vector_get(v,m++);
    }
  }

  int nind=int((par->ai).size());
  ln=0;
  for(int n=0;n<nind;n++){
    double p=pan2(nsnp,i0,L,(par->ai)[n],h1,J1);
    ln += -log(p);
  }
  ln /= nind;

  for(int l=0;l<L;l++)
    ln+= lambda*h1[l]*h1[l]/2;
  for(int i=0; i<nsnp; i++){
    if(i==i0) continue;
    for(int l=0;l<L*L;l++)
      ln+=lambda*J1[i][l]*J1[i][l]/2;
  }
  return ln;
}


void dlnl_psl(const gsl_vector *v,void *params,gsl_vector *df){   // first derivatives

  Param *par=(Param *)params;
  int L=par->L;
  int nsnp=(par->ai)[0].size();
  double lambda=par->lambda;
  int i0=par->i0;

  vector<double> s1(L);
  vector<vector<double> > s2(nsnp);

  vector<double> h1(L);
  vector<vector<double> > J1(nsnp);
  for(int i=0; i<nsnp; i++) J1[i].resize(L*L);

  int m=0;
  for(int l0=0;l0<L;l0++){
    h1[l0]=gsl_vector_get(v,m++);
    for(int i=0; i<nsnp; i++){ 
      if(i==i0) continue;
      for(int l1=0; l1<L; l1++)
        J1[i][L*l0+l1]=gsl_vector_get(v,m++);
    }
  }

  int nind=int((par->ai).size());

  for(int l=0;l<L;l++) s1[l]=0;
  for(int i=0; i<nsnp; i++){
    s2[i].resize(L*L);
    for(int l=0;l<L*L;l++) s2[i][l]=0;
  }

  for(int k=0;k<nind;k++){
    vector<double> peff(L);
    pan3(peff, nsnp, i0, L, (par->ai)[k], h1, J1);
    for(int l0=0;l0<L;l0++){
      double f=peff[l0]/nind;
      s1[l0]+=f;
      for(int j=0;j<nsnp;j++){
        if(j==i0) continue;
        short a=(par->ai)[k][j];
        if(a==0) continue;
        s2[j][L*l0+a-1]+=f;
      }
    }
  }

  for(int l0=0;l0<L;l0++){
    s1[l0] += lambda*h1[l0] - (par->f1)[l0];
    for(int j=0;j<nsnp;j++){
      if(j==i0) continue;
      for(int l1=0;l1<L;l1++)
        s2[j][L*l0+l1]+=-(par->f2)[j][L*l0+l1]
           +lambda*J1[j][L*l0+l1];
    }
  }

  m=0;
  for(int l0=0;l0<L;l0++){
    gsl_vector_set(df,m++,s1[l0]);
    for(int i=0; i<nsnp; i++){ 
      if(i==i0) continue;
      for(int l1=0;l1<L;l1++)
        gsl_vector_set(df,m++,s2[i][L*l0+l1]);
    }
  }
}

void ln_dln_psl(const gsl_vector *x,void *params,double *f,gsl_vector *df){

  *f=lnl_psl(x,params);
  dlnl_psl(x,params,df);

}

double lpr_psl(int i0, const vector<vector<short> > &ai, int L, double lambda, 
               vector<double> &h, vector<vector<double> > &J, int nprint, 
               unsigned int Imax, double Tol, int verbose){

  size_t iter=0;
  int status;

  int nind=int(ai.size());
  int nsnp=int(ai[0].size());
  vector<double> f1(nsnp);
  vector<vector<double> > f2(nsnp);

  f12(i0, ai, f1, f2, L);

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  int ndim = L+(nsnp-1)*L*L;

  my_func.n=ndim;         
  my_func.f=lnl_psl;
  my_func.df=dlnl_psl;
  my_func.fdf=ln_dln_psl;

  x=gsl_vector_alloc(ndim);
  T=gsl_multimin_fdfminimizer_vector_bfgs2;  // BFGS2 optimizer
  s=gsl_multimin_fdfminimizer_alloc(T,ndim);

  Param par={i0, ai, L, lambda, f1, f2};

  my_func.params=&par;
  gsl_vector_set_zero(x);  // initial guess

  gsl_multimin_fdfminimizer_set(s,&my_func,x,0.1,0.1);

  iter=0;
  do{
      iter++;
      status=gsl_multimin_fdfminimizer_iterate(s);
      if(iter%nprint==0 & verbose>0)
        cout << "  iteration # " << iter << ": " << s->f << endl;
      if(status){
        cerr << " GSL status code " << status << endl;
        exit(1);
      }
      status=gsl_multimin_test_gradient(s->gradient,Tol);
  }while(status==GSL_CONTINUE && iter< Imax); 
  if(iter==Imax)
    cerr << "BFGS2 iteration failed to converge after " 
         << Imax << " iterations\n";
  if(verbose > 0) cout << " Site " << i0+1 << ": " << iter 
       << " iterations, likelihood = " << s->f << endl;

  h.resize(L);
  J.resize(nsnp);
  double min=0;
  for(int i=0; i<nsnp; i++) J[i].resize(L*L);
  int m=0;
  for(int l0=0;l0<L;l0++){
    h[l0]=gsl_vector_get(s->x,m++);
    for(int i=0; i<nsnp; i++) for(int l1=0;l1<L;l1++){
        if(i==i0) J[i][L*l0+l1]=0;
        else
          J[i][L*l0+l1]=gsl_vector_get(s->x,m++);
    }
  }
  min=-nind*(s->f);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return min;
}


void f12(int i0, const vector<vector<short> > &si, vector<double> &f1,
         vector<vector<double> > &f2, int L){

  int n = si.size();
  int m = si[0].size();
  f1.resize(L);
  f2.resize(m);

  for(int i=0; i<m; i++){
    f2[i].resize(L*L);
    for(int l=0; l<L; l++) f1[l]=0; 
    for(int l0=0; l0<L; l0++) for(int l1=0; l1<L; l1++)
      f2[i][L*l0+l1]=0;
  }
  for(int k=0; k<n; k++){
    short a=si[k][i0];
    if(a==0) continue;
    f1[a-1]++;
    for(int j=0; j<m; j++){
      if(j==i0) continue;
      short b=si[k][j];
      if(b==0) continue;
      f2[j][L*(a-1)+b-1]++;
    }
  }
  for(int l0=0; l0<L; l0++){ 
    f1[l0]/=n;
    for(int j=0; j<m; j++){
      if(i0==j){
        for(int l1=0; l1<L; l1++){
          if(l0==l1)
            f2[j][L*l0+l1]=f1[l0];
          else
            f2[j][L*l0+l1]=0;
        }
      }
      else{
        for(int l1=0; l1<L; l1++)
          f2[j][L*l0+l1]/=n;
      }
    }
  }
}
