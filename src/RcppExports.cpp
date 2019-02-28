// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fillMyRSStable
NumericMatrix fillMyRSStable(NumericMatrix myRSStable, List& RSStriang, IntegerVector& myIndex, int m, float h);
RcppExport SEXP _strucchange_fillMyRSStable(SEXP myRSStableSEXP, SEXP RSStriangSEXP, SEXP myIndexSEXP, SEXP mSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type myRSStable(myRSStableSEXP);
    Rcpp::traits::input_parameter< List& >::type RSStriang(RSStriangSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type myIndex(myIndexSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< float >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(fillMyRSStable(myRSStable, RSStriang, myIndex, m, h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_strucchange_fillMyRSStable", (DL_FUNC) &_strucchange_fillMyRSStable, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_strucchange(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}