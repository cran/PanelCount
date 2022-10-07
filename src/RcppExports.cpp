// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// matVecProd
NumericMatrix matVecProd(NumericMatrix m, NumericVector v);
RcppExport SEXP _PanelCount_matVecProd(SEXP mSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(matVecProd(m, v));
    return rcpp_result_gen;
END_RCPP
}
// matVecProdSum
NumericMatrix matVecProdSum(NumericMatrix m, NumericVector ext, NumericVector v, NumericVector group);
RcppExport SEXP _PanelCount_matVecProdSum(SEXP mSEXP, SEXP extSEXP, SEXP vSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ext(extSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(matVecProdSum(m, ext, v, group));
    return rcpp_result_gen;
END_RCPP
}
// matVecProdSumExt
NumericMatrix matVecProdSumExt(NumericMatrix m, NumericVector ext, NumericVector ext2, NumericVector v, NumericVector group);
RcppExport SEXP _PanelCount_matVecProdSumExt(SEXP mSEXP, SEXP extSEXP, SEXP ext2SEXP, SEXP vSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ext(extSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ext2(ext2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(matVecProdSumExt(m, ext, ext2, v, group));
    return rcpp_result_gen;
END_RCPP
}
// groupProd
NumericMatrix groupProd(NumericVector v, NumericVector group);
RcppExport SEXP _PanelCount_groupProd(SEXP vSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(groupProd(v, group));
    return rcpp_result_gen;
END_RCPP
}
// groupSum
NumericMatrix groupSum(NumericVector v, NumericVector group);
RcppExport SEXP _PanelCount_groupSum(SEXP vSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(groupSum(v, group));
    return rcpp_result_gen;
END_RCPP
}
// groupSumMat
NumericMatrix groupSumMat(NumericMatrix m, NumericVector group);
RcppExport SEXP _PanelCount_groupSumMat(SEXP mSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(groupSumMat(m, group));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PanelCount_matVecProd", (DL_FUNC) &_PanelCount_matVecProd, 2},
    {"_PanelCount_matVecProdSum", (DL_FUNC) &_PanelCount_matVecProdSum, 4},
    {"_PanelCount_matVecProdSumExt", (DL_FUNC) &_PanelCount_matVecProdSumExt, 5},
    {"_PanelCount_groupProd", (DL_FUNC) &_PanelCount_groupProd, 2},
    {"_PanelCount_groupSum", (DL_FUNC) &_PanelCount_groupSum, 2},
    {"_PanelCount_groupSumMat", (DL_FUNC) &_PanelCount_groupSumMat, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_PanelCount(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
