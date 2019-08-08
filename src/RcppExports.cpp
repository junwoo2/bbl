// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pseudo_mle
List pseudo_mle(NumericMatrix xi, IntegerVector L, NumericVector Lambda, IntegerVector Nprint, IntegerVector Itmax, NumericVector Tol, IntegerVector Verbose);
RcppExport SEXP _BBM_pseudo_mle(SEXP xiSEXP, SEXP LSEXP, SEXP LambdaSEXP, SEXP NprintSEXP, SEXP ItmaxSEXP, SEXP TolSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Nprint(NprintSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Itmax(ItmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Tol(TolSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pseudo_mle(xi, L, Lambda, Nprint, Itmax, Tol, Verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BBM_pseudo_mle", (DL_FUNC) &_BBM_pseudo_mle, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_BBM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
