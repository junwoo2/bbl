const int Imax=10000;
const double Tol=1e-5;
struct Param{  
  int i0;   // central site index
  const std::vector<std::vector<short> > &ai;
  int L;
  int Lp;
  double lambda;
  const std::vector<double> &f1;
  const std::vector<std::vector<double> > &f2;
  double &lzp;
};

void f12(int i0, const std::vector<std::vector<short> > &si, 
         std::vector<double> &f1, std::vector<std::vector<double> > &f2, int L);

void pan3(std::vector<double> &peff, int nsnp, int i0, int L, int Lp,
          const std::vector<short> &ci, std::vector<double> h1, 
          const std::vector<std::vector<double> > &J1, double &lzp);

double lnl_psl(const gsl_vector *v, void *params);

void dlnl_psl(const gsl_vector *v, void *params, gsl_vector *df);

void ln_dln_psl(const gsl_vector *x, void *params, double *f, gsl_vector *df);

double lpr_psl(int i0, const std::vector<std::vector<short> > &si, int L, 
    double lambda, std::vector<double> &h, 
    std::vector<std::vector<double> > &J, int nprint,
    unsigned int imax, double tol, int verbose, double &lzp, bool numeric);

double pan2(int nsnp, int i0, int L, int Lp, const std::vector<short> &ci, 
    const std::vector<double> &h1, const std::vector<std::vector<double> > &J1,
    double &lzp);
