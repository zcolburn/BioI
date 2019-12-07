// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// euclidean_linker_cpp
Rcpp::NumericVector euclidean_linker_cpp(Rcpp::NumericMatrix input, double critDist, bool use_prog_bar);
RcppExport SEXP _Bioi_euclidean_linker_cpp(SEXP inputSEXP, SEXP critDistSEXP, SEXP use_prog_barSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type input(inputSEXP);
    Rcpp::traits::input_parameter< double >::type critDist(critDistSEXP);
    Rcpp::traits::input_parameter< bool >::type use_prog_bar(use_prog_barSEXP);
    rcpp_result_gen = Rcpp::wrap(euclidean_linker_cpp(input, critDist, use_prog_bar));
    return rcpp_result_gen;
END_RCPP
}
// find_min_dists_cpp
Rcpp::DataFrame find_min_dists_cpp(NumericMatrix mOne, NumericMatrix mTwo);
RcppExport SEXP _Bioi_find_min_dists_cpp(SEXP mOneSEXP, SEXP mTwoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mOne(mOneSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mTwo(mTwoSEXP);
    rcpp_result_gen = Rcpp::wrap(find_min_dists_cpp(mOne, mTwo));
    return rcpp_result_gen;
END_RCPP
}
